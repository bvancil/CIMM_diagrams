/*
  Routines for drawing CIMM dot diagrams

  A DotGroup can contain other DotGroups or DotClusters
*/
import stats; // for Gaussian random distribution
import graph; // for high-quality Circles and Ellipses;
// Global constants:
real phi=(1+sqrt(5))/2;
path unitCircle = Circle((0,0),1);

// utility functions:
path regular_star(int n, int jump) { // n points and jumping jump vertices (Note: jump=1 is a regular polygon.)
  pair[] z=new pair[n];
  for(int i=0; i<n; ++i) {
    z[i] = dir(90+360/n*i); // unit vector starting in up direction
  }
  z.cyclic = true;
  path p = z[0];
  for(int i=0; i<n; ++i) {
    p = p--z[i*jump];
  }
  p = p--cycle;
  return p;
}
path Ellipse(pair center, real a, real b) {
  return shift(center)*scale(a,b)*unitCircle;
}

// shapes:
// TODO: Allow arbitrary characters
int ball = 0; // lower-case oh
int square = 1; // square, U+25FB
int triangle_up = 2; // upward-pointing triangle
int triangle_down = 3; // downward-pointing triangle
int asterisk = 4; // 6-pointed cross
int star = 5; // a 5-pointed star
int pentagon = 6; // a regular pentagon
int hexagon = 7; // a regular hexagon
int heptagon = 8; // a regular heptagon
int octagon = 9; // a regular octagon
path[] shape_paths = {
  /* ball */ unitCircle,
  /* square */ polygon(4),
  /* triangle_up */ polygon(3),
  /* triangle_down */ invert*polygon(3),
  /* asterisk */ cross(6),
  /* star */ regular_star(5,2),
  /* pentagon */ polygon(5),
  /* hexagon */ polygon(6),
  /* heptagon */ polygon(7),
  /* octagon */ polygon(8)
};
path scaled_shape_path(int shape_name) { // give uniform horizontal size to shapes
    real desired_shape_size = 12pt;
    path shape = shape_paths[shape_name];
    pair shape_min = min(shape), shape_max = max(shape);
    real raw_horizontal_size = (shape_max-shape_min).x;
    return scale(desired_shape_size/raw_horizontal_size)*shape;
}
path shape_path(int shape_name) { // wrapper around scaled_shape_path in case no such shape exists
  if ((shape_name >= 0) && (shape_name < shape_paths.length)) return scaled_shape_path(shape_name);
  return scaled_shape_path(ball); // return default
}

struct DotCluster { // Homogeneous dots without a grouping boundary
  picture pic;
  pair center;
  real n;
  real a, b;  // horizontal and vertical, respectively, semiaxis lengths for the bounding ellipse
  int shape;
  real shape_size = 12pt; // TODO: fix to adjust depending on the shape
  pen base_p;
  pen dot_p;
  pair[] locations;

  bool _is_inside(pair p) { // Note that p is a relative location of a point.  The function will tell if it's inside the ellipse given by parameters a,b.
    // Note that we're contracting a & b by half the shape size so that it's more inside.
    return (length(scale(1/(this.a-this.shape_size/2),1/(this.b-this.shape_size/2))*p)<1);
  }
  bool _is_not_too_close(pair p) { // Note that p is a relative location of a point.  The function will tell if the location is not too close to previous locations.
    real tolerance = 1pt; // extra margin
    for(pair q : this.locations) {
      if (length(p-q)<this.shape_size+tolerance) return false;
    }
    return true;
  }
  real _primitive_random() { // return a random value between -1 and 1
    // Use a uniform distribution
    // return (2*unitrand()-1);
    // Use a clipped gaussian distribution to get more centrality
    real r;
    do {
      r = Gaussrand();
    } while (abs(r)>1);
    return r;
  }
  pair _random_location(real a=this.a, real b=this.b) {
    // Return a random location whose only constraint is on the range of each coordinate
    return ((a-this.shape_size/2)*this._primitive_random(), (b-this.shape_size/2)*this._primitive_random());
  }
  pair _semigood_random_location(real a=this.a, real b=this.b) {
    // Return a random location constrained by being not too close
    pair p;
    int i=0; // for debugging, count the iterations and warn
    do {
      p=this._random_location(a=a,b=b);
      assert(++i<100000, "taking too long to find a semigood random location"); // for debugging
    } while (!(this._is_not_too_close(p)));
    return p;
  }    
  pair _good_random_location(real a=this.a, real b=this.b) { // constrain by being inside and not too close
    pair p;
    int i=0; // for debugging, count the iterations and warn
    do {
      p=this._random_location(a=a,b=b);
      assert(++i<100000, "taking too long to find a good random location"); // for debugging
    } while (!(this._is_inside(p) && this._is_not_too_close(p)));
    return p;
  }
  pair _centroid() {
    pair sum;
    for (pair p : this.locations) sum += p;
    return sum/this.locations.length;
  }
  void _add_location() {
    // generate a decent point
    // add it to the list
    // TODO: grow by accretion
    // First, check if points in between will work.
    // Then, if they don't, step away from an edge point in a random direction.

    pair p;
    // For now, just generate it randomly inside the ellipse w/ a,b parameters
    //p = this._good_random_location();

    // Accretion code:
    if (this.locations.length == 0) p = (0,0); // if no points yet, start with zero;
    else {
      int i=0; // for debugging, count the iterations and warn
      pair c = this._centroid(); // middle of dots
      real dt = .02; // step size for accretion line of sight
      do { // repeat as many times as necessary to get a good trajectory that ends up inside the boundary ellipse
	pair q = this._semigood_random_location(a=5*this.a,b=5*this.b); // provisional point not too close
	//dot(q, blue+linewidth(0.5)); // for debugging
	// use line-of-sight accretion:
	for (real t=0; t<=1; t+=dt) {
	  p = (1-t)*q + t*c; // the point t along the segment from q to centroid c
	  //dot(dic,p, purple+linewidth(0.5)); // for debugging
	  if (!(this._is_not_too_close(p))) {
	    t -= dt; // back up
	    p = (1-t)*q + t*c; // recalculate last good p
	    break; // go with that p
	  }
	}
	assert(++i<10000, "taking too long to find an accretion trajectory");
      } while (!(this._is_inside(p)));
    }
    this.locations.push(p); // Add this decent point.
    //label(this.pic,format("%d",this.locations.length),p); // for debugging
  }
  void _set_pen() { // set the dot pen according to whether the dots are + or -.
    if (this.n > 0) {
      this.dot_p = this.base_p;
    } else if (this.n < 0) {
      this.dot_p = this.base_p + linetype(new real[] {0,1.5}); // fine dots
    }

  }
  picture draw(pen base_p=black) {
    // TODO: Dont' overwrite base_p if it has already been set.  
    this.base_p = base_p;
    this._set_pen();
    if (this.n<0) {
      for(pair p : this.locations) draw(this.pic,shift(p)*shape_path(this.shape), this.dot_p);
    }
    if (this.n>0) {
      for(pair p : this.locations) fill(this.pic,shift(p)*shape_path(this.shape), this.dot_p);
    }
    //draw(this.pic,box((-this.a,-this.b), (this.a,this.b)), green+white+dashed); // for debugging, draw random boundary box
    //draw(this.pic,Ellipse(this.center, this.a, this.b), red+white); // for debugging, draw ellipse for boundary
    //add(this.pic.fit());
    return this.pic;
  }
  void operator init(int n=0, int shape=ball, pair center=(0,0), pen base_p=black) {
    this.n = n;
    this.shape = shape;
    this.center = center;
    // Choose the pen:
    this.base_p = base_p;
    // For now I'll assume that the area will be proportional to the dots enclosed and confine dots there:
    real r = 10*sqrt(abs(n)); // length scale; TODO: adjust by scaling factor depending on the shape size
    // I'll also use an aspect ratio of the golden ratio phi (defined above):
    this.a = r*phi;
    this.b = r;
    srand(seconds()); // seed the random number generator
    for(int i=0; i<abs(this.n); ++i) { // choose where the balls go
      this._add_location();
    }
  }
}
from DotCluster unravel DotCluster;
void draw(DotCluster c, pen base_p=black) {
  // TODO: Dont' overwrite base_p if it has already been set.
  c.draw(base_p=base_p);
}

struct DotGroup {
  picture pic;
  pair new_center; // center for newly added children
  pair center; // center of this object
  real a, b;  // horizontal and vertical, respectively, semiaxis lengths for the bounding ellipse
  DotGroup[] child_groups;
  DotCluster[] child_clusters;
  bool boundary;
  string name;
  pen dot_p;
  pen base_p;
  pair _top(picture pic) { // return the top center coordinates of a picture
    pair s = size(pic);
    pair M = max(pic);
    return (M.x-s.x/2,M.y);
  }
  pair _bottom(picture pic) { // return the bottom center coordinates of a picture
    pair s = size(pic);
    pair m = min(pic);
    return (m.x+s.x/2,m.y);
  }
  pair _center(picture pic) { // return the center coordinates of a picture
    return (min(pic)+max(pic))/2;
  }
  void add_group(int n=1, DotGroup g) {
    this.child_groups.push(g);
    // Increment the size of the bounding ellipse in a heavyhanded way:
    // The 1.4 gives extra space for ellipse and label
    this.a += g.a*1.4; 
    this.b += g.b*1.4;
  }
  void add_dots(int n=1, int shape=ball, pair center=this.new_center, pen base_p=black) {
    DotCluster c = DotCluster(n=n, shape=shape, center=center, base_p=base_p);
    this.child_clusters.push(c);
    // Increment the size of the bounding ellipse in a heavyhanded way:
    this.a += c.a;
    this.b += c.b;
    this.new_center += scale(this.a,this.b)*dir(rand());
  }
  picture draw(pen base_p=black) {
    // TODO: Dont' overwrite base_p if it has already been set.
    // Draw contents:
    picture pic, pic_g, pic_c;
    for(DotGroup g : this.child_groups) add(pic_g,g.draw().fit(),this._bottom(pic_g),S);
    for(DotCluster c : this.child_clusters) add(pic_c,c.draw().fit(),this._bottom(pic_c),S);
    add(pic, pic_g.fit(), (0,0),E);
    add(pic, pic_c.fit(), (0,0),W);
    add(this.pic,pic.fit(),-this._center(pic));
    // Draw boundary:
    if (this.a==0) { // nothing inside this group; set a minimum size:
      this.b = 15pt;
      this.a = phi*this.b;
    }
    draw(this.pic,Ellipse(this.center,this.a,this.b));
    // Draw name/label:
    label(this.pic,this.name,shift(this.a,this.b)*this.center);
    //add(this.pic.fit());
    return this.pic;
  }
  void operator init(string name="", int n=0, int shape=ball, bool boundary=true, pen base_p=black) {
    this.name = name;
    this.boundary = boundary;
    // TODO: Are the following two lines unnecessary?
    this.a = 0; // set bounding ellipse to zero size before adding dots
    this.b = 0; // set bounding ellipse to zero size before adding dots
    this.add_dots(n=n, shape=shape, base_p=base_p);
    this.base_p = base_p;
  }
}
from DotGroup unravel DotGroup;
void draw(DotGroup g, pen base_p=black) {
  // TODO: Dont' overwrite base_p if it has already been set.
  g.draw(base_p=base_p);
}

// Usage: A DotGroup can contain dots or other DotGroups
// DotGroup U = DotGroup();
// U.add_group(n=3,DotGroup(n=2,shape=ball));
// DotGroup h = DotGroup(name="h");
// h.add_dots(n=4,shape=square);
// U.add_group(n=1,h);
// draw(U);

// Testing:

// Draw a single dot:
// DotCluster c = DotCluster(n=-10);
// draw(c,cmyk(blue));

// TODO: Fix overlap with squares and other shapes that aren't circles
// TODO: Fix filling of asterisk to filldraw of asterisk
// DotCluster c = DotCluster(n=-10,shape=square);
// draw(c,cmyk(blue));

// Draw a single homogeneous group of dots:
// DotGroup U = DotGroup(n=-10);
// draw(U,cmyk(blue));

// Draw a single inhomogeneous group of dots:
// DotGroup U = DotGroup(name="$U$", n=-10);
// U.add_dots(n=3);
// picture pic=U.draw(cmyk(blue));
// add(pic.fit());

// Draw a group of two homogeneous groups of dots:
DotGroup g = DotGroup(name="$g$", n=2, base_p=cmyk(blue));
DotGroup h = DotGroup(name="$h$", n=-3, base_p=cmyk(red));
DotGroup U = DotGroup(name="$U$", base_p=cmyk(black));
U.add_group(g);
U.add_group(g);
U.add_group(h);
add(U.draw().fit());
