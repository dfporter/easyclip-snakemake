import collections, pandas, os, re, glob, sys, importlib, pickle, subprocess, time, dill
from typing import List, Tuple, Union
from pprint import pprint
from pathlib import Path

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import scripts
import scripts.scheme
import scripts.makeBedDataFile
import scripts.scheme_signal_RNAs
import scripts.bedgraphs
import scripts.mapping

importlib.reload(scripts.scheme_signal_RNAs)
importlib.reload(scripts.makeBedDataFile)
importlib.reload(scripts.scheme)
importlib.reload(scripts.bedgraphs)
importlib.reload(scripts.mapping)

class exp(scripts.mapping.mappingMethods):
    """
    This is a holder of scheme, bed, count, and rna data locations.
    It differs from scheme_signal_RNAs objects, and the like, in that it
    only points to signal data locations, and never contains signal
    data, except during the creation of bed and scheme_signal_RNAs .data files, after
    which the large variables should go out of scope and be cleared from memory.
    
    This is intended to be used for dataset combinations using metaExp data objects.
    """
    
    def __init__(self, name: str = '0', file_paths: dict = {}):
        self.name = name
        self.file_paths = file_paths

    def read_scheme(self, fname: Union[str, None] = None) -> scripts.scheme.scheme:
        
        if fname is not None:
            self.file_paths['samples'] = fname
        self.scheme = scripts.scheme.scheme(self.file_paths['samples'])

        return self.scheme

    def make_signal_data_file(self, clobber=True, serialize=True):
        """ Make data/bed_x.data object of sequencing coverage over chromosome locations.
        Will serialize regardless of serialize argument."""
        
        data = self.file_paths['data']
        os.makedirs(data, exist_ok=True)

        if serialize:
            output_data_filename = f"{data}/bed_{self.name}.data"
            if (not clobber) and os.path.exists(output_data_filename):
                return       

        bedMaker = scripts.makeBedDataFile.bedDataFileMaker()
        
        # Originally, this function call was to dump to the .data file and returns a large object,
        # which we discarded. HTSeq.GenomicArray objects that are not empty will not serialize
        # currently on MacOS Mojave.
        """This fails with a seg fault on Mac in a Jupyter notebook (python3.8):
        _ga = HTSeq.GenomicArray('auto', stranded=True, typecode='i')
        _ga[HTSeq.GenomicInterval('1', 0, 100, '+')] += 1
        with open('./data/bed_m0.data', 'wb') as f:
            dill.dump(_ga, f)
        """
        print("Making bedgraphs objects from the folder {}, which has {} plus strand files.".format(
            self.file_paths['bedgraph'], len(glob.glob(self.file_paths['bedgraph'] + '/*.+.wig')))
             )
        
        # Returns the serializable form.
        self.beds = bedMaker.make_bed_data_file(
#            bed_folder=self.file_paths['beds'],
            bedgraphs_folder=self.file_paths['bedgraph'],
            output_data_filename=f'{data}/bed_{self.name}.data',
            use='read start',
            collapse_reads=False)
        
        return self.beds

    def make_scheme_signal_RNA_data_files(
        self,
        rna_data_object: Union[None, scripts.set_of_named_mRNAs.set_of_named_mRNAs] = None,
        no_clobber: bool = False, serialize: bool = True,
        counts_fname=None):
        
        if counts_fname is None:
            counts_fname = self.file_paths['top_dir'] + '/outs/counts/counts.bedgraphs.txt'
                          
        data = self.file_paths['data']
                          
        if rna_data_object is None:
            print(f"Loading RNA data file {data}/rna.data'...", end='')
            rna_data_object = pickle.load(open(f'{data}/rna.data', 'rb'))
            #rna_data_object.find_introns()
            print("Finished loading RNA data file.")

        if no_clobber and os.path.exists(f'{data}/signal_rnas_{self.name}.data'):
            print(f"scheme_signal_rna data file exists: {data}/signal_rnas_{self.name}.data." + \
                  " Not overwriting because no_clobber is on. ")
            if not(os.path.exists(counts_fname)):
                print(f"counts file {counts_fname} does not exist. Loading datafile to write it.")
                ssr = pickle.load(open(f'{data}/signal_rnas_{self.name}.data', 'rb'))
            else:  # Nothing to do, all outputs exist and no_clobber is on.
                return
        else:
            # Always returns the serializable form.
            if no_clobber and os.path.exists(f'{data}/bed_{self.name}.data'):
                # If serialize:
                print(f"make_scheme_signal_RNA_data_files(): Loading bed data file {data}/bed_{self.name}.data")
                self.beds = dill.load(open(f'{data}/bed_{self.name}.data', 'rb'))
                print("Finished loading bed data file.")
            else:
                self.beds = self.make_signal_data_file(clobber=not(no_clobber), serialize=serialize)

            print(f"make_scheme_signal_RNA_data_files(): Regenerating objects from serializable data.")
            self.beds.recover_bedgraph_objects()
            print("Finished regenerating objects.")

            ssr = scripts.scheme_signal_RNAs.scheme_signal_RNAs(
                set_of_bedgraphs=self.beds, set_of_RNAs=rna_data_object)

            ssr.add_scheme(self.file_paths['samples'])
            print(f"make_scheme_signal_RNA_data_files(): Making genomic array.")
            ssr.make_genomic_array(no_clobber=no_clobber)
            print(f"make_scheme_signal_RNA_data_files(): Assigning reads to genes...")
            ssr.assign_to_genes(verbose=True, no_clobber=no_clobber)
            print(f"Finished assigning reads to genes...")
                          
        #if 'counts' not in self.file_paths:
        #    self.file_paths['counts'] = self.file_paths['top_dir'] + '/outs/counts'

        if no_clobber or os.path.exists(counts_fname):
            print(f"Wrote raw counts file to {counts_fname}.")
            ssr.write_counts_file(
                fname=counts_fname, report_introns_separately=False)
        else:
            print(f"Did not write {counts_fname} because file exists and no_clobber is on.")
                          
        if serialize:
            if not(os.path.exists(f'{data}/signal_rnas_{self.name}.data')) or not(no_clobber):
                with open(f'{data}/signal_rnas_{self.name}.data', 'wb') as f:
                    pickle.dump(ssr, file=f)
    
    def duplication_rate(self) -> str:
        
        if not hasattr(self, 'scheme'):
            self.read_scheme()
        #beds = pickle.load(open('./data/bed_{0}.data'.format(self.name), 'rb'))
        
        #beds.duplication_rate()
        
        stats = scripts.bedAndFastaStats.bedAndFastaStats(
            scheme_object=self.scheme)

        reports = stats.duplication_rate(self.file_paths['beds'])
        
        if type(reports) != type(''):
            return 'No report\n {}'.format(reports)
        
        return reports
    


