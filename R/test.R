
defaults <- list(a=1, c=34)

for (arg in commandArgs(TRUE))
  eval(parse(text=arg))
  
for (nm in names(defaults))
  assign(nm, mget(nm, ifnotfound=list(defaults[[nm]]))[[1]])

print(a*2)
print(b*3)
print(c)

