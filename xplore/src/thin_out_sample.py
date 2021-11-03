# Take a nexus file containing trees and create a new file with every mth tree of the input file -- for thinning out long tree files from BEAST run where every tree is saved in file
__author__ = 'Lena Collienne'


def thin_out_sample(input_file, output_file, m):
    # Read sample from input_file and write into output_file every mth tree
    count = 0
    o_file = open(output_file, "w")
    with open(input_file) as i_file:
        for line in i_file:
            if count == 0:
                # Copy header + first tree from input to output file
                o_file.write(line)
                if re.match("tree", line):
                    # First line containing tree -- update counter
                    count += 1
            else:
                count +=1
                if (count-1) % m == 0:
                    o_file.write(line)
    o_file.write("End;")
    o_file.close()


if __name__ == '__main__':
    thin_out_sample("../simulations/posterior/primates/primates.trees", "../simulations/posterior/primates/primates_small.trees", 10000)
