#usage: python convert_cfg.py input.cfg > output.cfg

tokens = ('COLON','SEMICOLON','TERMINAL','NONTERMINAL')

t_COLON = r'\:'
t_SEMICOLON = r'\;'

t_NONTERMINAL = r'<[a-zA-Z_][a-zA-Z0-9_]*>'

def t_TERMINAL(t):
  r'\".*?\"'
  try:
    t.value = t.value[1:-1] #strip the quotes
  except ValueError:
    print("ValError: %s", t.value)
    t.value = ""
  return t

t_ignore = " \t\n"

def t_error(t):
  print("Illegal character '%s'" % t.value[0])
  t.lexer.skip(1)

import ply.lex as lex
lex.lex()


#from calclex import tokens

import ply.yacc as yacc

def p_cfg_rec(p) :
  'cfg : vardecl cfg'
  p[0] = p[1] + "\n" + p[2]

def p_cfg_end(p) :
  'cfg : vardecl'
  p[0] = p[1]

def p_empty(p):
    'empty :'
    pass

def p_vardecl(p) :
  'vardecl : NONTERMINAL rules'
  #p[0] = "<"+p[1]+">" + "::=" + p[2]
  p[0] = p[1] + "::=" + p[2]

def p_rules_rec(p) :
  'rules : rule rules'
  p[0] = p[1] + " | " + p[2]

def p_rules_end(p) :
  'rules : rule'
  p[0] = p[1] + ";"

def p_rule_empty(p) :
  'rule : COLON SEMICOLON'
#  p[0] = "\"<>\""
  p[0] = "\"()\""

def p_rule(p) :
  'rule : COLON symbollist SEMICOLON'
  p[0] = p[2]

def p_symbollist_rec(p) :
  'symbollist : symbol symbollist'
  p[0] = p[1] + " " + p[2]

def p_symbollist_end(p) :
  'symbollist : symbol'
  p[0] = p[1]

def p_symbol_NT(p) :
  'symbol : NONTERMINAL'
#  p[0] = "<"+p[1]+">"
  p[0] = p[1]

def p_symbol_Terminal(p) :
  'symbol : TERMINAL'
#  p[0] = "\"<"+p[1]+":1>\""
  p[0] = "\""+p[1]+"\""
def p_error(p):
    print "Syntax error in input!"

yacc.yacc()

import sys
try:
  f = open(sys.argv[1])
  p = yacc.parse(f.read())
  print p
except EOFError:
  print "EOF Err.."
