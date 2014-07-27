import itertools

persons = ['A','B','C']
products = ['p','q','r']

#trust and recommend relations with weights

#trust = {(0,1) : 0.3,
#         (0,2) : 0.7,
#         (1,0) : 0.2,
#         (1,2) : 0.6,
#         (1,3) : 0.2,
#         (4,2) : 0.2,
#         (3,4) : 0.5}

#rec = {(3,2) : 0.7,
#       (3,0) : 0.8,
#       (2,1) : 0.3,
#       (2,2) : 0.7}

# peer-to-peer trust values (<= 1)
trust = {(0,1) : 0.5,
         (0,2) : 0.2,
         (1,0) : 0.7,
         (2,0) : 0.6}

# sdsi-rep example:
#trust = {(0,1) : 0.7,
#         (1,2) : 0.9,
#         (2,0) : 0.5}


# direct recommendation weights (<=1)
rec = {(0,1) : 0.7,
       (1,1) : 0.8,
       (2,2) : 0.3}



def gen_system_free() :
  #for all pairs i,j in rang

  #generate equation for all pairs of persons e.g. <Trust_A_B> ::= "AB" | <Trust_A_C> <Trust_C_B> | ...
  for i,j in itertools.product(range(len(persons)),repeat=2) :
    joinpairs = ["<Trust_" + persons[i] + "_" + persons[x] + ">" + " <Trust_" + persons[x] + "_" + persons[j] + "> " for x in range(len(persons))]
    body = ""
    if len(joinpairs) > 0 :
      body = reduce(lambda x,y: x + " | " + y, joinpairs if trust.has_key((i,j)) else joinpairs[1:], ('"' + persons[i]+ "-" + persons[j] + '"') if trust.has_key((i,j)) else joinpairs[0])
    else :
      body = ('"' + persons[i]+ "-" + persons[j] + '"') if trust.has_key((i,j)) else ""
    
    if body == "" :
      body = "0";

    eq = "<Trust_" + persons[i] + "_" + persons[j] + "> ::= " + body + ";"
    print eq

  for i,j in itertools.product(range(len(persons)),range(len(products))) :
    joinpairs = ["<Trust_" + persons[i] + "_" + persons[x] + ">" + " <Rec_" + persons[x] + "_" + products[j] + "> " for x in range(len(persons))]
    body = ""
    if len(joinpairs) > 0 :
      body = reduce(lambda x,y: x + " | " + y, joinpairs if rec.has_key((i,j)) else joinpairs[1:], ('"' + persons[i]+"-" +products[j] + '"') if rec.has_key((i,j)) else joinpairs[0])
    else :
      body = ('"' + persons[i]+ "-" +products[j] + '"') if rec.has_key((i,j)) else ""
    if body == "" :
      body = "0";

    eq = "<Rec_" + persons[i] + "_" + products[j] + "> ::= " + body + ";"
    print eq


def gen_system_weighted() :
  #generate equation for all pairs of persons e.g. <Trust_A_B> ::= "AB" | <Trust_A_C> <Trust_C_B> | ...
  for i,j in itertools.product(range(len(persons)),repeat=2) :
    sum_of_p2p_weights = sum([trust[x,y] for x,y in trust.keys() if x==i])
    rest_weight = 1.0 - sum_of_p2p_weights
    if rest_weight > 0 :
      joinpairs = ['"' + str(rest_weight) + '"' + "<Trust_" + persons[i] + "_" + persons[x] + ">" + " <Trust_" + persons[x] + "_" + persons[j] + "> " for x in range(len(persons))]
    else :
      joinpairs = []
#    print str(i) + str(j) + " : " + str(joinpairs)
    body = ""
    if len(joinpairs) > 0 :
      if trust.has_key((i,j)) :
        body = reduce(lambda x,y: x + " | " + y, joinpairs, '"' + str(trust[i,j]) + '"')
        if trust[i,j] == 1.0:
          body = ('"' + str(trust[i,j]) + '"')
      else :
        body = reduce(lambda x,y: x + " | " + y, joinpairs[1:], joinpairs[0])
      #print body
    else :
      body = ('"' + str(trust[i,j]) + '"') if trust.has_key((i,j)) else ""
    if body == "" :
      body = "0";
    eq = "<Trust_" + persons[i] + "_" + persons[j] + "> ::= " + body + ";"
    print eq

  for i,j in itertools.product(range(len(persons)),range(len(products))) :
    sum_of_p2rec_weights = sum([rec[x,y] for x,y in rec.keys() if x==i])
    rest_rec = 1.0 - sum_of_p2rec_weights
    if rest_rec > 0 :
      joinpairs = ['"' + str(rest_rec) + '"' + "<Trust_" + persons[i] + "_" + persons[x] + ">" + " <Rec_" + persons[x] + "_" + products[j] + "> " for x in range(len(persons))]
    else :
      joinpairs = []
    body = ""
    if len(joinpairs) > 0 :
      if rec.has_key((i,j)) :
        body = reduce(lambda x,y: x + " | " + y, joinpairs, '"' + str(rec[i,j]) + '"')
      else :
        body = reduce(lambda x,y: x + " | " + y, joinpairs[1:], joinpairs[0])
    else :
      body = ('"' + str(rec[i,j]) + '"') if rec.has_key((i,j)) else ""
    if body == "" :
      body = "0";
    eq = "<Rec_" + persons[i] + "_" + products[j] + "> ::= " + body + ";"
    print eq



if __name__ == "__main__":
  #import sys
  #foo(sys.argv[1])
  #gen_system_free()
  gen_system_weighted()

