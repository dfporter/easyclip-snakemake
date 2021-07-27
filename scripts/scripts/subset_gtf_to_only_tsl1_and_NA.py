import re, sys


def subset_to_only_tsl1_and_NA(input_file, output_file):
    outf = open(output_file, 'w')

    with open(input_file) as f:
        for li in f:
            
            if (('transcript_support_level "1"' in li) or ('transcript_support_level "NA"' in li)):
                
                outf.write(li)

    outf.close()
    

if __name__ == '__main__':
    
    input_file, output_file = sys.argv[1], sys.argv[2]
   
    subset_to_only_tsl1_and_NA(input_file, output_file)