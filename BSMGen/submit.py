import os
import argparse

if __name__ == '__main__':

    parser = argparse.ArgumentParser()    
    parser.add_argument('-i', '--input', type=str, help="Input pythia card", required=True)
    parser.add_argument('-j', '--njob', type=int, default = 10, help="Number of jobs")
    parser.add_argument('-q', '--queue', type=str, default = "1nh", help="LSFBATCH queue name")
    parser.add_argument('--joblabel', type=str, default = "LCDgun", help="job label")

    args = parser.parse_args()

    # existing files
    os.system("ls /eos/cms/store/group/phys_susy/razor/tchannel/Delphes/%s | grep root  > %s_existing.list" %(args.joblabel,args.joblabel))
    f = open("%s_existing.list" %args.joblabel)
    existing = f.readlines()
    f.close()

    if not os.path.isdir(args.joblabel): 
        os.system("mkdir %s" %args.joblabel)
    if not os.path.isdir('/eos/cms/store/group/phys_susy/razor/tchannel/Delphes/'+args.joblabel):
        os.system("mkdir /eos/cms/store/group/phys_susy/razor/tchannel/Delphes/%s/\n" %args.joblabel)
    if not os.path.isdir('/eos/cms/store/group/phys_susy/razor/tchannel/ntuples/'+args.joblabel):
        os.system("mkdir /eos/cms/store/group/phys_susy/razor/tchannel/ntuples/%s/\n" %args.joblabel)
    if not os.path.isdir('/eos/cms/store/group/phys_susy/razor/tchannel/log/'+args.joblabel):
        os.system("mkdir /eos/cms/store/group/phys_susy/razor/tchannel/log/%s/\n" %args.joblabel)
    os.system("chmod 777 /eos/cms/store/group/phys_susy/razor/tchannel/Delphes/%s/\n" %args.joblabel) # So that batch can write, apparently
    for i in range(args.njob):
        myfile = "%s_%i.root\n" %(args.joblabel, i)
        if myfile in existing: continue
        script = open("%s/%s_%i.src" %(args.joblabel, args.joblabel, i), "w")
        script.write("cd /afs/cern.ch/user/q/qnguyen/BSMatLHC/BSMGen\n")
        script.write("eval `scramv1 run -sh`\n")
        #script.write("cd -\n")     
        script.write("cd %s\n" %os.getcwd())
        script.write("source setup.sh\n")
        script.write("mkdir /tmp/qnguyen/%s/\n" %args.joblabel)
        script.write("./GenPythia %s /tmp/qnguyen/%s/%s_%i \n" %(args.input, args.joblabel, args.joblabel, i))
        script.write("cd /afs/cern.ch/user/q/qnguyen/BSMatLHC/delphes\n")
        script.write("./DelphesHepMC cards/CMS_PhaseII/CMS_PhaseII_0PU.tcl /tmp/qnguyen/%s/%s_%i.root /tmp/qnguyen/%s/%s_%i.hepmc\n" %(args.joblabel, args.joblabel, i,args.joblabel, args.joblabel, i)) 
        script.write("cd /afs/cern.ch/user/q/qnguyen/BSMatLHC/BSMApp\n")
        script.write("cp /tmp/qnguyen/%s/%s_%i.root /eos/cms/store/group/phys_susy/razor/tchannel/Delphes/%s\n" %(args.joblabel, args.joblabel, i,args.joblabel))
        script.write("./CMSApp  /tmp/qnguyen/%s/%s_%i.root --output=/tmp/qnguyen/%s/%s_%i.root --tchannel --delphes\n" % (args.joblabel, args.joblabel, i, args.joblabel, args.joblabel, i))
        script.write("cp /tmp/qnguyen/%s/%s_%i.root /eos/cms/store/group/phys_susy/razor/tchannel/ntuples/%s\n" %(args.joblabel, args.joblabel, i,args.joblabel))
        script.write("rm /tmp/qnguyen/%s/%s_%i.hepmc\n" %(args.joblabel, args.joblabel, i))
        script.write("rm /tmp/qnguyen/%s/%s_%i.root\n" %(args.joblabel, args.joblabel, i))
        script.write("rm /tmp/qnguyen/%s/%s_%i.root\n" %(args.joblabel, args.joblabel, i))
        script.close()
        os.system("bsub -q %s -o /eos/cms/store/group/phys_susy/razor/tchannel/log/%s/%s_%i.log -J %s_%i < %s/%s_%i.src" %(args.queue, args.joblabel, args.joblabel, i, args.joblabel, i, args.joblabel, args.joblabel, i));
        print "Submitting job %i to the queue %s...\n" %(i,args.queue)
