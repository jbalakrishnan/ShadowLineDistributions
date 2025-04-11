load("shadow.sage")
import statistics
from collections import Counter


def parse_and_filter_dist(input_file, filter = True):
    """
    Reads in input_file containing all shadow line data (e.g., 1483.a1_3/1483.a1_3_out) and parses it in the following way:
    (1) If filter = True, produces an output file with slopes corresponding to fields K with R that is p-divisible or p dividing h_K omitted
    (i.e., whenever a line of input_file contains the string 'R is p-divisible!' or 'p divides h_K !', the next several lines through the end of that
    particular example are deleted).
    (2) Splits off the set of all slopes and writes it to an output file (separate from the one produced in (1), if that step is invoked).
    (3) If filter = True, returns the filtered slope distribution.
    """
    skipping = False
    if filter == True:
        # remove the entries with R that is p-divisible or p dividing h_K
        with open(input_file, 'r', encoding='utf-8') as infile, open(input_file+"_filtered_temp", 'w', encoding='utf-8') as outfile:
            target_string = ["R is p-divisible!", "p divides h_K !", "stats:  "]
            for line in infile:
                if skipping:
                    if "--------------------" in line:
                        skipping = False
                else:
                    if any(marker in line for marker in target_string):
                        skipping = True
                    else:
                        outfile.write(line)
        outfile.close()

    with open(input_file, 'r', encoding='utf-8') as infile, open(input_file+"_just_slopes", 'w', encoding='utf-8') as outfile:
         for line in infile:
            if "slope list: " in line:
                outfile.write(line)

    if filter == True:
        # now process the filtered slope file
        slopes = []
        target_substring = "slope:  "
        with open(input_file+"_filtered_temp", 'r', encoding='utf-8') as infile:
            for line in infile:
                if target_substring in line:
                    start_index = line.find(target_substring)
                    slope_value = line[start_index + len(target_substring):].strip()
                    slopes.append(slope_value)
        slope_new = [sage_eval(x) for x in slopes]
        p = slope_new[0].parent().prime()
        return partp(slope_to_modp(slope_new,p),p)


def slope_to_modpsquare_nonan(slopelist,p):
    """
    For a list of shadow line slopes with (E, p) nonanomalous, return the list of values mod p^2.
    """
    modpsquare = []
    Zn = Zmod(p^2)
    for x in slopelist:
        try:
            modpsquare.append(Zn(x))
        except ValueError:
            modpsquare.append(['oo', Zn(1/x)])
    return modpsquare


def partp_nonan(modlist, p):
    """
    Count each of the mod p^2 slope values mod p, and for each of the values mod p, refine this by
    producing the list giving the distribution of the next coefficient of the slope (mod p^2).

    For example, for (709.a1, 3), we have the following:
    59 slopes that are (0, 1) mod 3,
    50 slopes that are (1,1) mod 3,
    61 slopes that are (2, 1) mod 3, and
    55 slopes that are (1,0) mod 3,
    producing the  distribution [59, 50, 61, 55].

    This function produces the output
    [[18, 18, 23], [14, 17, 19], [19, 18, 24], [14, 18, 23]]
    where [18, 18, 23] gives the frequency of the 59 slopes that are (0, 1), (3, 1), and (6, 1) mod 9, respectively.
    """
    lists = [[modlist.count(i + p*a) for a in range(p)] for i in range(p)]
    modinfty = [a[1] for a in modlist if type(a) == list]
    listinfty = [modinfty.count(0 + p*a) for a in range(p)]
    return lists + [listinfty]


def slope_to_modpsquare_anomalous(slopelist, p, invert=False):
    """
    For a list of shadow line slopes with (E, p) anomalous, return the list of values mod p^2.
    """
    modpsquare = []
    Zn = Zmod(p^2)
    if invert==True:
        slopelist = [1/x for x in slopelist]
    for x in slopelist:
        try:
            modpsquare.append(ZZ(Zn(x)))
        except ValueError:
            modpsquare.append(['oo'])
    return modpsquare


def modpsquare_partp(modlist, mode, p):
    """
    For the modal value mod p, producing the list giving the distribution of the next coefficient of the slope (mod p^2).
    """
    restricted_modlist = []
    for x in modlist:
        if type(x) != list:
            if (x-mode) %p == 0:
                restricted_modlist.append((ZZ(x)-ZZ(mode))/ZZ(p))
    return [restricted_modlist.count(i) for i in range(p)]

def separating_n(slopelist):
    """
    Given a list of shadow line slopes mod p^N, returns the minimal n such that all slopes mod p^n are distinct.
    """
    p = slopelist[0].parent().prime()
    for prec in range(1,30):
        slopelisttrunc = [x + O(p^prec) for x in slopelist]
        k = Counter(slopelisttrunc).keys()
        v = Counter(slopelisttrunc).values()
        if list(v) == sum(v)*[1]:
            return prec
            break

print(75*"=")
print("Nonanomalous data")
print(75*"=")
#this is the list of nonanomalous distributions
for x in [['709.a1',3],['997.c1',3],['1627.a1',3],['2677.a1',3],['709.a1',5],['1531.a1',5],['1621.a1',5],['1873.a1',5],['1907.a1',5],['1933.a1',5],['643.a1',7],['709.a1',7],['997.c1',7],['1613.a1',7],['1627.a1',7]]:
    print("(E, p) = (%s, %s)"%(x[0],x[1]))
    y = x[0] + '_%s'%(x[1])
    l = '../../nonanomalous/' + y +  '/' + y + '_out' #change the path as appropriate
    parse_and_filter_dist(l)

    with open('../../nonanomalous/' + y +  '/' + y + '_out_just_slopes', 'r', encoding='utf-8') as infile: #change the path as appropriate
        for line in infile:
            if "slope list:  " in line:
                start_index = line.find("slope list:  ")
                slopes = line[len("slope list:  "):].strip()
        slopes = sage_eval(slopes)
    print("Separated slopes mod p^n with n = ", separating_n(slopes))
    modpsquare_slopes = slope_to_modpsquare_nonan(slopes,x[1])
    print("Mod p^2 distribution: ", partp_nonan(modpsquare_slopes,x[1]))
    print(50*'*')

print(75*"=")
print("Anomalous data")
print(75*"=")
#this is the list of anomalous distributions
for x in [['433.a1', 3], ['643.a1', 3],['1058.a1', 3],['1483.a1', 3],['1613.a1', 3],['1933.a1', 3],['6293.d1', 3],['36781.b1', 3],['433.a1', 5],['563.a1', 5],['997.c1', 5],['6011.a1', 7],['2251.a1', 11],['1933.a1', 13],['709.a1', 29],['1483.a1', 31]]:
    print("(E, p) = (%s, %s)"%(x[0],x[1]))
    y = x[0] + '_%s'%(x[1])
    l = '../../anomalous/' + y +  '/' + y + '_out'  #change the path as appropriate
    filtered = parse_and_filter_dist(l)


    #computing the mode of the unfiltered anomalous distribution
    with open('../../anomalous/' + y +  '/' + y + '_out_just_slopes', 'r', encoding='utf-8') as infile:  #change the path as appropriate
        for line in infile:
            if "slope list:  " in line:
                start_index = line.find("slope list:  ")
                slopes = line[len("slope list:  "):].strip()
        slopes = sage_eval(slopes)
        slopes_mod_p = slope_to_modp(slopes,x[1])
        slope_mode = statistics.mode(slopes_mod_p)
    print("Separated slopes mod p^n with n = ", separating_n(slopes))
    print("Distribution of slopes mod p: ", partp(slopes_mod_p,x[1]))
    print("Distribution of filtered slopes mod p: ", filtered)
    if (x[0] == '433.a1' and x[1] == 3) or (x[0] == '6293.d1' and x[1] == 3):
        inv = True
        slope_mode = 0
        print("Will work with inverse slopes since mode is +Infinity")
    else:
        inv = False
    print("Slopes mod p mode: ", slope_mode)
    if x[0] not in ['1058.a1', '6293.d1']:
        modpsquare_slopes = slope_to_modpsquare_anomalous(slopes,x[1],inv)
        print("Mod p^2 distribution: ", modpsquare_partp(modpsquare_slopes,slope_mode,x[1]))
    print(50*'*')


"""
Here is the output that is produced by the above:

===========================================================================
Nonanomalous data
===========================================================================
(E, p) = (709.a1, 3)
Separated slopes mod p^n with n =  10
Mod p^2 distribution:  [[18, 18, 23], [14, 17, 19], [19, 18, 24], [14, 18, 23]]
**************************************************
(E, p) = (997.c1, 3)
Separated slopes mod p^n with n =  10
Mod p^2 distribution:  [[22, 19, 29], [24, 20, 20], [23, 22, 20], [21, 15, 21]]
**************************************************
(E, p) = (1627.a1, 3)
Separated slopes mod p^n with n =  10
Mod p^2 distribution:  [[21, 24, 24], [23, 20, 11], [26, 21, 17], [19, 14, 21]]
**************************************************
(E, p) = (2677.a1, 3)
Separated slopes mod p^n with n =  11
Mod p^2 distribution:  [[16, 15, 19], [15, 24, 23], [14, 20, 18], [24, 18, 19]]
**************************************************
(E, p) = (709.a1, 5)
Separated slopes mod p^n with n =  7
Mod p^2 distribution:  [[4, 10, 12, 11, 11], [11, 12, 11, 9, 9], [6, 12, 12, 4, 4], [6, 7, 6, 14, 11], [8, 8, 9, 12, 6], [8, 7, 12, 9, 5]]
**************************************************
(E, p) = (1531.a1, 5)
Separated slopes mod p^n with n =  8
Mod p^2 distribution:  [[9, 2, 6, 7, 12], [4, 6, 8, 12, 12], [10, 9, 9, 6, 12], [5, 8, 8, 16, 7], [9, 11, 11, 7, 6], [3, 9, 9, 8, 13]]
**************************************************
(E, p) = (1621.a1, 5)
Separated slopes mod p^n with n =  7
Mod p^2 distribution:  [[8, 12, 13, 5, 5], [7, 6, 10, 11, 5], [15, 12, 11, 8, 11], [8, 9, 12, 10, 8], [9, 10, 8, 13, 9], [9, 9, 8, 4, 9]]
**************************************************
(E, p) = (1873.a1, 5)
Separated slopes mod p^n with n =  7
Mod p^2 distribution:  [[16, 13, 8, 12, 10], [11, 10, 7, 11, 4], [11, 6, 7, 10, 9], [9, 7, 17, 8, 9], [10, 6, 8, 12, 9], [7, 8, 7, 7, 8]]
**************************************************
(E, p) = (1907.a1, 5)
Separated slopes mod p^n with n =  8
Mod p^2 distribution:  [[8, 7, 10, 12, 6], [3, 1, 12, 9, 9], [3, 8, 10, 11, 7], [9, 9, 5, 2, 7], [4, 4, 6, 12, 8], [7, 13, 8, 6, 9]]
**************************************************
(E, p) = (1933.a1, 5)
Separated slopes mod p^n with n =  6
Mod p^2 distribution:  [[7, 7, 7, 9, 9], [7, 12, 8, 9, 11], [5, 12, 5, 7, 7], [9, 11, 8, 12, 8], [10, 14, 10, 9, 14], [15, 7, 10, 16, 7]]
**************************************************
(E, p) = (643.a1, 7)
Separated slopes mod p^n with n =  5
Mod p^2 distribution:  [[2, 4, 4, 3, 3, 2, 6], [2, 7, 8, 5, 3, 3, 3], [3, 2, 1, 2, 3, 6, 7], [6, 2, 6, 8, 2, 3, 2], [4, 6, 7, 5, 2, 7, 3], [7, 4, 8, 6, 4, 1, 4], [2, 6, 9, 2, 5, 4, 5], [1, 5, 4, 7, 5, 1, 3]]
**************************************************
(E, p) = (709.a1, 7)
Separated slopes mod p^n with n =  5
Mod p^2 distribution:  [[3, 1, 7, 3, 4, 0, 6], [6, 6, 6, 3, 4, 5, 3], [5, 5, 8, 3, 9, 5, 5], [5, 6, 3, 3, 2, 6, 3], [6, 2, 3, 2, 6, 7, 3], [2, 8, 2, 3, 3, 3, 3], [4, 5, 3, 4, 8, 5, 4], [3, 4, 3, 7, 4, 7, 5]]
**************************************************
(E, p) = (997.c1, 7)
Separated slopes mod p^n with n =  6
Mod p^2 distribution:  [[4, 3, 7, 8, 4, 3, 4], [4, 3, 4, 3, 4, 7, 2], [4, 2, 3, 4, 1, 5, 5], [8, 3, 6, 5, 5, 3, 7], [5, 7, 5, 2, 4, 4, 4], [5, 3, 5, 3, 1, 2, 3], [3, 9, 4, 3, 3, 2, 5], [4, 5, 3, 1, 5, 3, 3]]
**************************************************
(E, p) = (1613.a1, 7)
Separated slopes mod p^n with n =  5
Mod p^2 distribution:  [[2, 5, 8, 7, 5, 10, 7], [7, 6, 4, 8, 3, 7, 6], [8, 2, 5, 10, 6, 5, 7], [2, 2, 2, 2, 3, 6, 6], [4, 5, 4, 8, 6, 3, 3], [0, 2, 5, 4, 5, 4, 5], [10, 5, 4, 3, 1, 4, 5], [6, 2, 4, 6, 10, 6, 7]]
**************************************************
(E, p) = (1627.a1, 7)
Separated slopes mod p^n with n =  6
Mod p^2 distribution:  [[3, 8, 6, 6, 2, 4, 5], [8, 3, 8, 4, 6, 4, 8], [7, 6, 7, 2, 5, 8, 4], [4, 6, 12, 5, 8, 4, 8], [4, 6, 2, 4, 3, 4, 3], [4, 6, 4, 5, 4, 5, 5], [6, 5, 4, 8, 8, 7, 6], [8, 2, 4, 4, 5, 2, 5]]
**************************************************
===========================================================================
Anomalous data
===========================================================================
(E, p) = (433.a1, 3)
Separated slopes mod p^n with n =  8
Distribution of slopes mod p:  [25, 21, 26, 208]
Distribution of filtered slopes mod p:  [14, 18, 19, 126]
Will work with inverse slopes since mode is +Infinity
Slopes mod p mode:  0
Mod p^2 distribution:  [71, 60, 77]
**************************************************
(E, p) = (643.a1, 3)
Separated slopes mod p^n with n =  11
Distribution of slopes mod p:  [25, 28, 139, 36]
Distribution of filtered slopes mod p:  [0, 0, 102, 0]
Slopes mod p mode:  2
Mod p^2 distribution:  [42, 47, 50]
**************************************************
(E, p) = (1058.a1, 3)
Separated slopes mod p^n with n =  7
Distribution of slopes mod p:  [23, 25, 20, 25]
Distribution of filtered slopes mod p:  [0, 7, 0, 0]
Slopes mod p mode:  +Infinity
**************************************************
(E, p) = (1483.a1, 3)
Separated slopes mod p^n with n =  10
Distribution of slopes mod p:  [32, 147, 28, 29]
Distribution of filtered slopes mod p:  [0, 110, 0, 0]
Slopes mod p mode:  1
Mod p^2 distribution:  [36, 61, 50]
**************************************************
(E, p) = (1613.a1, 3)
Separated slopes mod p^n with n =  10
Distribution of slopes mod p:  [24, 164, 31, 50]
Distribution of filtered slopes mod p:  [0, 133, 0, 0]
Slopes mod p mode:  1
Mod p^2 distribution:  [45, 62, 57]
**************************************************
(E, p) = (1933.a1, 3)
Separated slopes mod p^n with n =  10
Distribution of slopes mod p:  [43, 24, 170, 33]
Distribution of filtered slopes mod p:  [0, 0, 125, 0]
Slopes mod p mode:  2
Mod p^2 distribution:  [59, 57, 54]
**************************************************
(E, p) = (6293.d1, 3)
Separated slopes mod p^n with n =  7
Distribution of slopes mod p:  [23, 21, 22, 46]
Distribution of filtered slopes mod p:  [12, 15, 17, 0]
Will work with inverse slopes since mode is +Infinity
Slopes mod p mode:  0
**************************************************
(E, p) = (36781.b1, 3)
Separated slopes mod p^n with n =  15
Distribution of slopes mod p:  [33, 24, 116, 19]
Distribution of filtered slopes mod p:  [0, 0, 84, 0]
Slopes mod p mode:  2
Mod p^2 distribution:  [40, 34, 42]
**************************************************
(E, p) = (433.a1, 5)
Separated slopes mod p^n with n =  7
Distribution of slopes mod p:  [21, 8, 13, 193, 11, 16]
Distribution of filtered slopes mod p:  [0, 0, 0, 175, 0, 0]
Slopes mod p mode:  3
Mod p^2 distribution:  [38, 37, 36, 46, 36]
**************************************************
(E, p) = (563.a1, 5)
Separated slopes mod p^n with n =  8
Distribution of slopes mod p:  [14, 17, 170, 16, 10, 8]
Distribution of filtered slopes mod p:  [0, 0, 151, 0, 0, 0]
Slopes mod p mode:  2
Mod p^2 distribution:  [32, 34, 34, 37, 33]
**************************************************
(E, p) = (997.c1, 5)
Separated slopes mod p^n with n =  8
Distribution of slopes mod p:  [10, 17, 23, 15, 192, 14]
Distribution of filtered slopes mod p:  [0, 0, 0, 0, 166, 0]
Slopes mod p mode:  4
Mod p^2 distribution:  [50, 36, 35, 33, 38]
**************************************************
(E, p) = (6011.a1, 7)
Separated slopes mod p^n with n =  6
Distribution of slopes mod p:  [13, 9, 11, 7, 226, 8, 10, 5]
Distribution of filtered slopes mod p:  [0, 0, 0, 0, 213, 0, 0, 0]
Slopes mod p mode:  4
Mod p^2 distribution:  [24, 37, 28, 37, 41, 27, 32]
**************************************************
(E, p) = (2251.a1, 11)
Separated slopes mod p^n with n =  5
Distribution of slopes mod p:  [2, 1, 3, 3, 2, 2, 2, 3, 181, 2, 4, 0]
Distribution of filtered slopes mod p:  [0, 0, 0, 0, 0, 0, 0, 0, 179, 0, 0, 0]
Slopes mod p mode:  8
Mod p^2 distribution:  [19, 15, 14, 16, 17, 17, 11, 14, 22, 18, 18]
**************************************************
(E, p) = (1933.a1, 13)
Separated slopes mod p^n with n =  5
Distribution of slopes mod p:  [2, 4, 2, 2, 1, 1, 1, 4, 3, 2, 1, 8, 229, 4]
Distribution of filtered slopes mod p:  [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 222, 0]
Slopes mod p mode:  12
Mod p^2 distribution:  [12, 26, 20, 15, 23, 9, 21, 18, 18, 17, 18, 22, 10]
**************************************************
(E, p) = (709.a1, 29)
Separated slopes mod p^n with n =  4
Distribution of slopes mod p:  [0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 196, 0, 0, 2, 0, 0]
Distribution of filtered slopes mod p:  [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 196, 0, 0, 0, 0, 0]
Slopes mod p mode:  24
Mod p^2 distribution:  [3, 13, 9, 8, 7, 7, 5, 7, 6, 6, 4, 2, 7, 7, 5, 8, 6, 8, 7, 7, 7, 9, 7, 3, 4, 7, 7, 12, 8]
**************************************************
(E, p) = (1483.a1, 31)
Separated slopes mod p^n with n =  4
Distribution of slopes mod p:  [1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 196, 0, 0, 1, 0, 0, 2, 0, 0, 0, 0, 0]
Distribution of filtered slopes mod p:  [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 195, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
Slopes mod p mode:  20
Mod p^2 distribution:  [7, 7, 5, 8, 7, 3, 4, 3, 7, 8, 10, 5, 8, 7, 11, 4, 9, 7, 3, 10, 6, 9, 3, 8, 4, 4, 3, 8, 7, 7, 4]
**************************************************
"""


