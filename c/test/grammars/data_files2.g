  <datafile>    ::= <record> . [ <datafile> ]
  <record>      ::= <field> [ ; <record> ]
  <field>       ::= <integer>  |  <string>
  <integer>     ::= <digit> [ <integer> ]
  <string>      ::= " <stringchars> "
  <stringchars> ::= <char> [ <stringchars> ]
  <char>        ::= <letter>  |  <digit>
<letter>      ::= A | B | C | E | G | H | J | K | L | M | N | 
                           P | R | S | T | V | W | X | Y | Z
<digit>        ::= 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9

