B ::= 1 | B VarDeclaration
VarDeclaration ::= Type Identifier ";"
Type ::= "int" "[" "]" | "boolean" | "int" | Identifier
Identifier ::= "<IDENTIFIER>"
