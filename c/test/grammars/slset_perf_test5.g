<x0> ::= "<b:1>" <x0><x1> | "<a:2>" <x1><x2> | "<>";
<x1> ::= "<a:2>" <x2> | "<a:5,b:1>" <x4> | "<>";
<x2> ::= "<a:2,b:1>" <x2><x3> | "<b:3>" <x2><x2> | "<a:2>";
<x3> ::= "<a:3,b:1>" <x2><x2> | "<a:2,b:2>" <x3><x3> | "<>";
<x4> ::= "<a:1,b:2>" <x2><x5> | "<a:2,b:3>" <x1><x4> | "<>";
<x5> ::= "<a:3,b:2>" <x2><x2> | "<a:1,b:3>" <x1><x1> | "<>";
<x6> ::= "<a:1,b:3>" <x6><x2> | "<a:4,b:2>" <x6><x7> | "<>";
<x7> ::= "<a:2,b:5>" <x7><x3> | "<a:13,b:6>" <x5><x7> | "<>";
