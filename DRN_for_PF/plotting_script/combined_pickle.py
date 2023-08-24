import os, sys
import pickle

file1 = sys.argv[1]
file2 = sys.argv[2]

db = {}
  
myfile = open(file1,"rb")
db[os.path.splitext(file1)[0]]= pickle.load(myfile)
myfile.close()
print(file1)

myfile = open(file2,"rb")
db[os.path.splitext(file2)[0]]= pickle.load(myfile)
myfile.close()
print(file2)

print(db)
myfile = open("tmp.pkl","wb")
pickle.dump(db, myfile)
myfile.close()
