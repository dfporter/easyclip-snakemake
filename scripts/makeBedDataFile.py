import glob, importlib, pickle, os, dill
import scripts
import scripts.bedgraphs
import scripts.collapse_duplicates
import subprocess
from subprocess import PIPE

importlib.reload(scripts.bedgraphs)
importlib.reload(scripts.collapse_duplicates)


def out(command):
    #result = subprocess.run(command, stdout=PIPE, stderr=PIPE, universal_newlines=True, shell=True)
    #return result.stdout
    print(command)
    os.system(command)


def _mk(dirname):
    if not os.path.exists(dirname):
        os.system('mkdir ' + dirname)

class bedDataFileMaker():
    
    def __init__(self):
        pass
    
    def make_bed_data_file(
        self, bedgraph_list=None, bed_folder=None,
        bedgraphs_folder=None,
        use='read start', verbose=True, output_data_filename='data/bed_181005.data',
        collapse_reads=False):
        
        #bed_dir = '/Users/dfporter/pma/miseq/Runs/171210/sams/consensus/'
        #bedgraph_list = glob.glob(
        #    '/Users/dfporter/pma/miseq/Runs/171210/sams/consensus/' + '/{prefix}*_deletions.+.wig'.format(prefix='ATCGTG')
        #    )  # PCBP1 in 171210
        if collapse_reads:
            print('Making backups of bed files...')
            _mk('backup/')
            _mk('backup/bed/')
            _mk('backup/bed/{0}/'.format(os.path.basename(bed_folder.rstrip('/'))))
            out('rsync -av {0} {1}'.format(bed_folder, os.path.basename(bed_folder.rstrip('/'))))
            print('...Made backups. Collapsing...')
           #print('Initial folder size: ')
            #out('du -h {0}'.format(bed_folder))

            scripts.collapse_duplicates.collapse_bed_dir(bed_folder, 'temp_bed/')
            os.system('mv temp_bed/*bed ' + bed_folder + '/')
            #print('Folder size after collapsing:')
            #print(out('du -h {0}'.format(bed_folder)))
        
        if bedgraphs_folder is not None and (bed_folder is not None):
            raise ValueError("Cannot set both bedgraphs_folder and bed_folder")
        
        if bed_folder is not None:
            beds = scripts.bedgraphs.set_of_bedgraphs(
                bed_folder=bed_folder,
                use='read start',
                verbose=True,
            #    bedgraphs_folder=bed_dir
                )
        if bedgraphs_folder is not None:
            beds = scripts.bedgraphs.set_of_bedgraphs(
                bedgraphs_folder=bedgraphs_folder,
                verbose=True
            )
        #print(beds.__dict__)
        
        print(f"Made {len(beds.bedgraphs)} bedgraph data objects.")
        
        if not os.path.exists('data/'):
            os.system('mkdir data/')
        
        print("bedDataFileMaker(): Attempting to serialize the data... This may take a long time...")
        beds.make_serializable()
        
        print(f"bedDataFileMaker(): Writing to {output_data_filename}... This may take a long time...")
        with open(output_data_filename, 'wb') as fh:
            dill.dump(beds, file=fh)
        fh.close()
        
        print(f"bedDataFileMaker(): Wrote bed/bedgraph data file to {output_data_filename}.")
        
        return beds