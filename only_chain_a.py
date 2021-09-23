
import glob

def main():
    pdb_file = glob.glob("*.pdb")
    num = 0
    for pdb in pdb_file:
        f = open(str(num) + ".pdb","w")
        num1 = 0
        for i in open(pdb):
            ff = i.split()
            if len(ff) > 5:
                if ff[4] == "B":
                    num1 += 1
                    i = i.replace(" B ", "  A ")
                    i = i.replace(f[1], str(TEST(f[1], num1)))
            else:
                f.write(i)
        f.close()
        num += 1

def TEST(pdb, num):
    num1 = num
    for i in range(len(str(pdb) - len(str(num)))):
        num1 = " " + num1
    return num1



if __name__ == '__main__':
    main()
