"""
I mean to codify a 1-D textual dot notation and write a parser to turn that into a picture, maybe through TikZ or maybe Asymptote.


Using http://en.wikipedia.org/wiki/Extended_Backus%E2%80%93Naur_Form (EBNF):

digit excluding zero = "1" | "2" | "3" | "4" | "5" | "6" | "7" | "8" | "9" ;
digit                = "0" | digit excluding zero ;
natural number = digit excluding zero, { digit }, natural number ;
integer = "0" | [ "-" ], natural number ;
symbol = "" | "" | "" | ""; (*FIXME! Really any unicode symbol that is not combining or a number or just shapes? *)
dot cluster = integer , symbol
dot cluster sequence = "" | dot cluster , [(" "|"+"),
group = "(", dot cluster sequence, ")"
signed group = ("" | "-" | integer), group
signed group sequence = 

"""
