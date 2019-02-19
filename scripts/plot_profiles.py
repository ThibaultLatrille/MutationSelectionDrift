#!python3
import argparse
import os


def plot_profiles(args_input, args_infer, args_output):
    nperline = 53

    for profile in args_infer:
        prefs = args_output + "/" + os.path.basename(profile).replace(".siteprofiles", ".prefs")
        with open(profile, 'r') as r:
            with open(prefs, 'w') as w:
                r.readline()
                w.write("# POSITION WT SITE_ENTROPY PI_A PI_C PI_D PI_E PI_F PI_G PI_H PI_I PI_K PI_L PI_M PI_N PI_P PI_Q PI_R PI_S PI_T PI_V PI_W PI_Y\n")
                for line in r:
                    line_split = line.replace('\n', '').split('\t')
                    w.write(line_split[0] + " A 1.0 " + " ".join(line_split[1:]) + "\n")

        os.system("dms_logoplot {0} {0}.pdf --nperline {1}".format(prefs, nperline))
        diff_file = prefs + ".diff"
        os.system("dms_merge {0} sum {1} --minus {2}".format(diff_file, prefs, args_input))
        os.system("dms_logoplot {0} {0}.pdf --nperline {1}".format(diff_file, nperline))
        os.system('dms_correlate {1} {0} {0}.corr --name1 "Initial" --name2 "Estimated"'.format(prefs, args_input))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--input', required=True, type=str, dest="input")
    parser.add_argument('-p', '--infer', required=True, type=str, nargs='+', dest="infer")
    parser.add_argument('-o', '--output', required=True, type=str, dest="output")
    args = parser.parse_args()
    plot_profiles(args.input, args.infer, args.output)
