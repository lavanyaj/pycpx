import fileinput



value = None
sym = None
enum = None
meaning = None

StatusCodes = {}
t = 0
for line in fileinput.input():
    if t == 0:
        value = int(line.rstrip())
        pass
    elif t == 1:
        sym = str(line.rstrip())
        pass
    elif t == 2:
        enum = str(line.rstrip())
        pass
    elif t == 3:
        meaning = str(line.rstrip())
        StatusCodes[value] = {'sym': sym, 'enum': enum, 'meaning':meaning}
        pass
    t = (t + 1)%4
    pass
print len(StatusCodes)
codeValues = sorted(StatusCodes.keys())

for v in codeValues:
    print "#define %s %d" % (StatusCodes[v]['sym'], v)
    pass

print "switch(status) {"
for v in codeValues:
    sym = StatusCodes[v]['sym']
    enum = StatusCodes[v]['enum']
    meaning = StatusCodes[v]['meaning']
    print "case %s:" % enum
    print "  return Status(\"%s\", %s);" % (meaning, sym)
    pass
print "default:"
print "  return Status(\"Unknown status code from Cplex\");"
print "}"


