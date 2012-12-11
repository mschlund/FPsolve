<datafile>    ::= <record> { <record> }
<record>      ::= <field> { ; <field> } .
<field>       ::= <integer>  |  <string>
<integer>     ::= <digit> { <digit> }
<string>      ::= " <char> { <char> } "
<char>        ::= <letter>  |  <digit>
<letter>      ::= A | B | C | E | G | H | J | K | L | M | N | 
                           P | R | S | T | V | W | X | Y | Z
<digit>       ::= 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9

