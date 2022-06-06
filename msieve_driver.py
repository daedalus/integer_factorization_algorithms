
def msieve_factor(n):
  import subprocess, re
  BIN = "/home/dclavijo/code/FACTORING/msieve/msieve"
  #ARGS="-v %d" % n
  tmp = []
  proc = subprocess.Popen([BIN,"-t","8","-v",str(n)],stdout=subprocess.PIPE)
  for line in proc.stdout:
    #the real code does filtering here
    #print("test:", line.rstrip())
    line = line.rstrip().decode("utf8")
    #print(line)
    if re.search("factor: ",line):
      tmp += [int(line.split()[2])]
  return tmp

#n = 196616705776852626592749272255227899
import sys
n = int(sys.argv[1])
print(msieve_factor(n))
