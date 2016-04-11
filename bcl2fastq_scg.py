#!/usr/bin/python

#!usr/bin/python

import re
import sys
import subprocess
import os
import datetime
import time
import argparse
import fnmatch
import shutil
from Bio.Seq import Seq

class FlowcellLane:
    
    def __init__(self, run_name, lane_index, bcl2fastq_version, lims_url, lims_token, year, month,
                 rev_complement=True):

        # For not just get/put everything in properties
        #self.lane_project_dxid = self.properties['lane_project_dxid']
        self.run_name = run_name
        self.lane_index = lane_index
        self.lims_url = lims_url
        self.lims_token = lims_token
        self.bcl2fastq_version = bcl2fastq_version
        self.year = year
        self.month = month
        self.rev_complement = rev_complement
        
        self.sample_sheet = None
        self.use_bases_mask = None
        self.flowcell_id = None
        #self.output_dir = None

        # Get flowcell id
        run_elements = self.run_name.split('_')
        flowcell_info = run_elements[3]
        self.flowcell_id = flowcell_info[1:6]

    def describe(self):
        print "Sequencing run: %s" % self.run_name
        print "Flowcell lane index: %s" % self.lane_index

    def create_sample_sheet(self):
        '''
        create_sample_sheet.py -r ${seq_run_name} \
            -t ${UHTS_LIMS_TOKEN} \
            -u ${UHTS_LIMS_URL} \
            -b 2 \
            -l ${SGE_TASK_ID}
        '''

        ## Create samplesheet
        sample_sheet = '%s_L%s_samplesheet.csv' % (self.run_name, self.lane_index)
        if os.path.exists(sample_sheet):
            self.sample_sheet = sample_sheet
            print 'Sample sheet %s already exists. Skipping creation.' % sample_sheet
            return sample_sheet
        else:
            program_path = '/srv/gsfs0/software/gbsc/scgpm_lims/current/scgpm_lims/scripts/create_sample_sheet.py'
            command = 'python %s -r %s -t %s -u %s -b %d -l %d' % (program_path, self.run_name, self.lims_token, self.lims_url, self.bcl2fastq_version, self.lane_index)
            stdout,stderr = self._createSubprocess(cmd=command, pipeStdout=True)
            self.sample_sheet = '%s_L%d_samplesheet.csv' % (self.run_name, self.lane_index)
            stdout_elements = stdout.split()
            self.sample_sheet = stdout_elements[1]
            print 'This is the sample_sheet: %s' % self.sample_sheet
        
    def get_use_bases_mask(self):
        '''
        command = "python calculate_use_bases_mask.py {runinfoFile} {sampleSheet} {lane}"
        gbsc_utils.createSubprocess(cmd=command)
        '''
        
        run_info_file = 'RunInfo.xml'
        program_path = '/srv/gsfs0/projects/gbsc/workspace/pbilling/code/calculate_use_bases_mask_2.py'
        command = 'python %s %s %s %s' % (program_path, run_info_file, self.sample_sheet, self.lane_index)
        stdout,stderr = self._createSubprocess(cmd=command, pipeStdout=True)
        self.use_bases_mask = stdout
        print 'This is use_bases_mask value: %s' % self.use_bases_mask

        use_bases_mask_file = 'use_bases_mask.txt'
        with open(use_bases_mask_file, 'w') as OUT:
            OUT.write(self.use_bases_mask)
        return self.use_bases_mask

    def run_bcl2fastq(self):
        '''
        bcl2fastq --output-dir ${new_run_dir}/${seq_run_name}/Unaligned_L${SGE_TASK_ID} \
            --sample-sheet ${new_run_dir}/${seq_run_name}/${seq_run_name}_L${SGE_TASK_ID}_samplesheet.csv \
            --ignore-missing-bcls \
            --ignore-missing-filter \
            --ignore-missing-positions \
            --barcode-mismatches 1 \
            --use-bases-mask ${SGE_TASK_ID}:Y*,n*,Y*
        '''

        self.output_dir = 'Unaligned_L%d' % self.lane_index
        
        command = 'bcl2fastq '
        command += '--output-dir %s ' % self.output_dir
        command += '--sample-sheet %s ' % self.sample_sheet
        command += '--barcode-mismatches %d ' % 1
        command += '--use-bases-mask %d:%s ' % (int(self.lane_index), self.use_bases_mask)
        command += '--ignore-missing-bcls '
        command += '--ignore-missing-filter '
        command += '--ignore-missing-positions'

        stdout,stderr = self._createSubprocess(cmd=command, pipeStdout=True)

    def publish_fastqs(self, run_dir, pub_run_dir):
        # Make published lane directory if does not already exist
        pub_lane_dir = pub_run_dir + '/' + 'L' + str(self.lane_index)
        if not os.path.exists(pub_lane_dir):
            os.makedirs(pub_lane_dir)

        # Change to current lane directory
        lane_dir = run_dir + '/Unaligned_L' + str(self.lane_index)
        flowcell_dir = lane_dir + '/' + self.flowcell_id
        reports_dir = lane_dir + '/' + 'Reports'

        # Rename and copy all fastq files from lane directory
        for file in os.listdir(lane_dir):
            if fnmatch.fnmatch(file, '*.fastq.gz'):
                new_file_name = self._rename_fastq(file)
                src = lane_dir + '/' + file
                dst = pub_lane_dir + '/' + new_file_name
                shutil.copy(src, dst)

        for file in os.listdir(flowcell_dir):
            if fnmatch.fnmatch(file, '*.fastq.gz'):
                new_file_name = self._rename_fastq(file)
                src = flowcell_dir + '/' + file
                dst = pub_lane_dir + '/' + new_file_name
                shutil.copy(src, dst)

        # Copies all files in the Reports folder
        src = reports_dir
        dst = pub_lane_dir + '/Reports'
        shutil.copytree(src, dst)

        #cp -r Unaligned_L${SGE_TASK_ID}/Reports/* \
        #    ${pub_run_dir}/${year}/${month}/${seq_run_name}/L${SGE_TASK_ID}

    def tar_files(self, pub_run_dir):
        pub_lane_dir = pub_run_dir + '/L' + str(self.lane_index)
        os.chdir(pub_lane_dir)

        tar_file = '%s_L%d_pf.fastq.tar' % (self.run_name, self.lane_index)
        cmd = 'tar -cvf %s *.fastq.gz' % tar_file
        self._createSubprocess(cmd)

    def make_html(self, pub_run_dir):

        result_html_filename = 'L%d.html' % self.lane_index
        result_html_path = '%s/L%d/%s' % (pub_run_dir, self.lane_index, result_html_filename)
        index_link = 'http://scg-data.stanford.edu/PublishedResults/%d/%s/%s/L%d/Reports/html/index.html' % (self.year, self.month, self.run_name, self.lane_index)
        tar_link = 'http://scg-data.stanford.edu/PublishedResults/%d/%s/%s/L%d/%s_L%d_pf.fastq.tar' % (self.year, self.month, self.run_name, self.lane_index, self.run_name, self.lane_index)
        tar_path = '%s/L%d/%s_L%d_pf.fastq.tar' % (pub_run_dir, self.lane_index, self.run_name, self.lane_index)

        html_code = '''<html>
    <head>
         <title>Results for %s Lane %d</title>
    </head>
    <body>
    <h2>Results for %s Lane %d</h2>
    <li><a href="%s">Demultiplexing Statistics</a></li>
    <table border="1" cellspacing="0" cellpadding="4">
        <tr>
            <td>Description</td>
            <td>Barcode</td>
            <td>Read Number</td>
            <td>Format</td>
            <td>File Size</td>
            <td>Link</td>
        </tr>
        <tr>
            <td>Unmapped Post-Filter Reads</td>
            <td>All</td>
            <td>Both</td>
            <td>FASTQ</td>
            <td align=right>NA</td>
            <td><a href="%s">Download</a></td>
        </tr>
    </body>
</html>''' % (self.run_name, self.lane_index, self.run_name, self.lane_index, index_link, tar_link)

        # QA check to make sure that tarball exists
        if os.path.isfile(tar_path):
            with open(result_html_path, 'w') as HTML:
                HTML.write(html_code)
        else:
            print 'Could not find tarball: %s' % tar_path
            sys.exit()

    def _rename_fastq(self, fastq_filename):

        elements = fastq_filename.split('_')
        if len(elements) < 5 or len(elements) > 7:
            print 'WARNING: fastq filename has unusual number of elements : %s' % fastq
            sys.exit()
        
        # Two barcodes : lane1_TCTCGCGC_TCAGAGCC_S47_L001_R2_001.fastq.gz
        elif len(elements) == 7:
            lane = elements[0]
            barcodes = [elements[1], elements[2]]
            read = elements[5]

            # Need to use reverse complement of second barcode
            if self.rev_complement == True:
                rev_complement = self._reverse_complement(barcodes[1])
                new_barcodes = [barcodes[0], rev_complement]
            elif self.rev_complement == False:
                new_barcodes = [barcodes[0], barcodes[1]]

            lane_index_match = re.match(r'lane(\d)', lane)
            if lane_index_match:
                lane_index = lane_index_match.group(1)
            else:
                print 'Could not determine lane index: %s' % lane
                sys.exit()
            read_index_match = re.match(r'R(\d)', read)
            if read_index_match:
                read_index = read_index_match.group(1)
            else:
                print 'Could not determine read index: %s' % read
            
            # new format : 151202_BRISCOE_0270_BC847TACXX_L1_TAGGCATG-GGCTCTGA_1_pf.fastq.gz
            new_fastq_file = '%s_L%s_%s-%s_%s_pf.fastq.gz' % (self.run_name, lane_index, new_barcodes[0], new_barcodes[1], read_index)  
        
        # One barcode : lane1_TCTCGCGC_S47_L001_R2_001.fastq.gz
        elif len(elements) == 6:
            lane = elements[0]
            barcode = elements[1]
            read = elements[4]

            lane_index_match = re.match(r'lane(\d)', lane)
            if lane_index_match:
                lane_index = lane_index_match.group(1)
            else:
                print 'Could not determine lane index: %s' % lane
                sys.exit()
            read_index_match = re.match(r'R(\d)', read)
            if read_index_match:
                read_index = read_index_match.group(1)
            else:
                print 'Could not determine read index: %s' % read
                sys.exit()

            # new format : 151202_BRISCOE_0270_BC847TACXX_L1_TAGGCATG_1_pf.fastq.gz
            new_fastq_file = '%s_L%s_%s_%s_pf.fastq.gz' % (self.run_name, lane_index, barcode, read_index)
        
        # No barcode : Undetermined_S1_L001_R1_001.fastq.gz
        elif len(elements) == 5:
            lane = elements[2]
            read = elements[3]

            lane_index_match = re.match(r'L00(\d)', lane)
            if lane_index_match:
                lane_index = lane_index_match.group(1)
            else:
                print 'Could not determine lane index: %s' % lane
                sys.exit()
            read_index_match = re.match(r'R(\d)', read)
            if read_index_match:
                read_index = read_index_match.group(1)
            else:
                print 'Could not determine read index: %s' % read
            
            # new format : 151106_LYNLEY_0515_AC7F31ACXX_L1_unmatched_1_pf.fastq.gz
            new_fastq_file = '%s_L%s_unmatched_%s_pf.fastq.gz' % (self.run_name, lane_index, read_index)
        else:
            print "Could not get metadata for:\nfastq: %s\nlane: %s\nrun: %s" % (fastq, lane, run)
            pass
        return new_fastq_file

    def _reverse_complement(self, sequence):
        dna = Seq(sequence)
        rev_complement = dna.reverse_complement()
        return rev_complement

    def _createSubprocess(self, cmd, pipeStdout=False, checkRetcode=True):
        """
        Function : Creates a subprocess via a call to subprocess.Popen with the argument 'shell=True', and pipes stdout and stderr. Stderr is always  piped, but stdout can be turned off.
                 If the argument checkRetcode is True, which it is by defualt, then for any non-zero return code, an Exception is
                             raised that will print out the the command, stdout, stderr, and the returncode when not caught. Otherwise, the Popen instance will be return, in which case the caller must
                           call the instance's communicate() method (and not it's wait() method!!) in order to get the return code to see if the command was a success. communicate() will return
                             a tuple containing (stdout, stderr). But at that point, you can then check the return code with Popen instance's 'returncode' attribute.
        Args     : cmd   - str. The command line for the subprocess wrapped in the subprocess.Popen instance. If given, will be printed to stdout when there is an error in the subprocess.
                             pipeStdout - bool. True means to pipe stdout of the subprocess.
                             checkRetcode - bool. See documentation in the description above for specifics.
        Returns  : A two-item tuple containing stdout and stderr, respectively.
        """
        stdout = None
        if pipeStdout:
            stdout = subprocess.PIPE
            stderr = subprocess.PIPE
        popen = subprocess.Popen(cmd,shell=True,stdout=stdout,stderr=subprocess.PIPE)
        if checkRetcode:
            stdout,stderr = popen.communicate()
            if not stdout: #will be None if not piped
                stdout = ""
            stdout = stdout.strip()
            stderr = stderr.strip()
            retcode = popen.returncode
            if retcode:
                #below, I'd like to raise a subprocess.SubprocessError, but that doens't exist until Python 3.3.
                raise Exception("subprocess command '{cmd}' failed with returncode '{returncode}'.\n\nstdout is: '{stdout}'.\n\nstderr is: '{stderr}'.".format(cmd=cmd,returncode=retcode,stdout=stdout,stderr=stderr))
            return stdout,stderr
        else:
            return popen
            #return stdout,stderr

def main():
    sys.path.append('/srv/gsfs0/software/gbsc/scgpm_lims/current/scgpm_lims/scripts')
    
    parser = argparse.ArgumentParser(description='Convert raw bcl files to fastq and demultiplex')
    parser.add_argument('--run_name', '-r', dest='run_name', required=True, type=str, help='string with name of sequencing run')
    parser.add_argument('--lane_index', '-l', dest='lane_index', required=True, type=int, help='int with year of sequencing run processing')
    parser.add_argument('--bcl2fastq_version', '-b', dest='bcl2fastq_version', required=True, type=int, help='string with month of sequencing run processing')
    parser.add_argument('--lims_url', '-u', dest='lims_url', required=True, type=str)
    parser.add_argument('--lims_token', '-t', dest='lims_token', required=True, type=str)
    parser.add_argument('--year', '-y', dest='year', required=True, type=int)
    parser.add_argument('--month', '-m', dest='month', required=True, type=str)
    parser.add_argument('--no_reverse_complement', '-n', dest='rev_complement', action='store_false', required=False)
    args = parser.parse_args()

    print 'Performing bcl2fastq conversion & demultiplexing with arguments:\n'
    print args

    if not len(sys.argv) > 1:
        parser.print_help()
        sys.exit()

    new_runs_dir = '/srv/gsfs0/projects/seq_center/Illumina/RunsInProgress'
    pub_runs_dir = '/srv/gsfs0/projects/seq_center/Illumina/PublishedResults'
    run_dir = new_runs_dir + '/' + args.run_name
    pub_run_dir = pub_runs_dir + '/' + str(args.year) + '/' + args.month + '/' + args.run_name
    os.chdir(run_dir)

    lane = FlowcellLane(args.run_name, args.lane_index, args.bcl2fastq_version, args.lims_url, args.lims_token, args.year, args.month, args.rev_complement)
    lane.describe()
    print 'Creating sample sheet\n'             # Tested
    sample_sheet = lane.create_sample_sheet()
    print 'Get use bases mask\n'                # Tested
    use_bases_mask = lane.get_use_bases_mask()
    print 'Convert bcl to fastq files'
    lane.run_bcl2fastq()                        # Not tested
    print 'Copy renamed fastq files to published runs directory'
    lane.publish_fastqs(run_dir, pub_run_dir)   # Tested
    print 'Tar renamed files in published directory for easy download'
    lane.tar_files(pub_run_dir)                 # Tested
    print 'Generate HTML file for lane, to download tar files'
    lane.make_html(pub_run_dir)                 # Tested

if __name__ == '__main__':
    main()
