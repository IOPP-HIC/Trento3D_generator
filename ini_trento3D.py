# generate minimum bias events and collect the entropy information in every task 
import sys
import h5py as h5
import pandas as pd
import multiprocessing as mp
from glob import glob
import numpy as np
import os
import math
import cmath
from subprocess import check_output, call
import time


def run_jet_position(comment,cpuid,threadid,projectile,target, num_of_event, reduced_thickness, beam_energy,cross_section,\
        bmin,bmax,grid_max,grid_step,eta_max,eta_step,\
        mean_coeff,std_coeff,skew_coeff,skew_type,jacobian,fluctuation,\
        nucleon_width,constituent_width,constituent_num,min_dist):
    proc=call(['./trento3d_subnucleon/build/src/trento3d', '%s'%projectile, '%s'%target, '%s'%num_of_event, "-q" ,   "-b", '-p', '%s'%reduced_thickness,\
         '-e', '%s'%beam_energy, '-x', '%s'%cross_section, '-o', './events/%s/%s_%s.h5'%(comment,cpuid,threadid), \
         '--b-min','%s'%bmin,'--b-max','%s'%bmax, '--xy-max','%s'%grid_max,'--xy-step','%s'%grid_step, '--eta-max','%s'%eta_max, '--eta-step', '%s'%eta_step,\
        '--mean-coeff','%s'%mean_coeff,'--std-coeff','%s'%std_coeff,'-t','%s'%skew_coeff,'-r','%s'%skew_type,'-j','%s'%jacobian,'-k','%s'%fluctuation,\
         '-w','%s'%nucleon_width, '-v', '%s'%constituent_width, '-m', '%s'%constituent_num, '-d','%s'%min_dist])
    
        
def get_multi_list_and_centrality_cut(path,dir_num,grid_max,grid_step,eta_max,eta_step,comment,cpuid):
    mult=[]
    
    
    ngridxy  = int(2 * grid_max / grid_step)
    ngrideta = int(2 * eta_max  / eta_step)+1
    if not os.path.exists('./events/%s/mult_seta_max_%s.dat'%(comment,cpuid)):
        all_events=h5.File('./events/%s/all_events.h5'%comment,'w')
        count=0
        for threadid in range(dir_num):
            flist=h5.File('./events/%s/%s_%s.h5'%(comment,cpuid,threadid),'r')
            print(threadid)
            for key in flist.keys():
                seta=flist[key]['matter_density'][...].reshape(ngridxy, ngridxy,ngrideta).sum(axis=0).sum(axis=0)
                temp=seta.max()
                mult.append( temp )

                all_events.create_dataset('event_%s'%count, data=flist[key]['matter_density'][...].flatten(), dtype="f", compression="gzip", compression_opts=5)
                count=count+1
            flist.close()
        np.savetxt('./events/%s/mult_seta_max_%s.dat'%(comment,cpuid),mult)
        all_events.close()

        #for threadid in range(dir_num):
        #    call('rm -rf ./events/%s/%s_%s.h5'%(comment,cpuid,threadid),shell=True)


    else:
        mult=np.loadtxt('./events/%s/mult_seta_max_%s.dat'%(comment,cpuid))

    mult_sort = mult.sort(reverse=True)
    
    def get_limit_centrality(mult_sort, cent=10):    
        nindex = int(len(mult_sort) * cent / 100)
        cent_low_limit = mult_sort[nindex]
        return cent_low_limit

    with open('./events/%s/centrality_cut_seta_max_%s.dat'%(comment,cpuid),"w") as f:
        for i in range(0, 100):
            cent_low_limit = get_limit_centrality(mult_sort,cent=i)
            f.write('{}% {}\n'.format(i, cent_low_limit))


def get_multi_list_and_centrality_cut2(path,dir_num,grid_max,grid_step,eta_max,eta_step,comment,cpuid):
    mult=[]
    
    
    ngridxy  = int(2 * grid_max / grid_step)
    ngrideta = int(2 * eta_max  / eta_step)+1

    if not os.path.exists('./events/%s/all_events.h5'%comment):
        outputpath = './events/%s/all_events.h5'%comment
        count=0
        for threadid in range(dir_num):
            inputpath = './events/%s/%s_%s.h5'%(comment,cpuid,threadid)
            flist=h5.File(inputpath,'r')
            keys_tem = list(flist.keys())
            flist.close()
            for key in keys_tem:

                sourcepath = "{0}/{1}".format(key,"matter_density")
                targetpath = 'event_%s'%count

                call("h5copy -i {0} -o {1} -s {2} -d {3}".format(inputpath,outputpath,sourcepath,targetpath),shell=True)
                print("h5copy -i {0} -o {1} -s {2} -d {3}".format(inputpath,outputpath,sourcepath,targetpath))
                count=count+1
    if not os.path.exists('./events/%s/mult_seta_max_%s.dat'%(comment,cpuid)):
        flist=h5.File('./events/%s/all_events.h5'%comment,'r')
        for ikey in range(len(flist.keys())):
            print (ikey)
            seta=flist["event_%d"%ikey][...].reshape(ngridxy, ngridxy,ngrideta).sum(axis=0).sum(axis=0)
            temp=seta.max()
            mult.append( temp )

        flist.close()
        np.savetxt('./events/%s/mult_seta_max_%s.dat'%(comment,cpuid),mult)
        mult = np.array(mult)
    else:
        mult=np.loadtxt('./events/%s/mult_seta_max_%s.dat'%(comment,cpuid))
        #mult = list(mult)

    print (mult)
    mult_sort = np.sort(mult)[::-1]
    
    def get_limit_centrality(mult_sort, cent=10):    
        nindex = int(len(mult_sort) * cent / 100)
        cent_low_limit = mult_sort[nindex]
        return cent_low_limit

    with open('./events/%s/centrality_cut_seta_max_%s.dat'%(comment,cpuid),"w") as f:
        for i in range(0, 100):
            cent_low_limit = get_limit_centrality(mult_sort,cent=i)
            f.write('{}% {}\n'.format(i, cent_low_limit))
    


   
def collect_event(centlist,comment,cpuid,nchoose=3000):

    def get_centlimit(fpath, cent='0_10'):
        c1, c2 = [i + '%' for i in cent.split('_')]
        with open(fpath) as f:
            d = dict(line.split() for line in f.readlines())
            d['100%'] = 0
        return (float(d[c1]), float(d[c2]))

    fpath='./events/%s/centrality_cut_seta_max_%s.dat'%(comment,cpuid)
    all_events_path='./events/%s/all_events.h5'%comment
    cent_events_path='./events/%s/cent_events.h5'%comment
    
    cent_event_file=h5.File(cent_events_path,'w')

    eventlist=np.array([],dtype=int)
    for cent in centlist:
        
        cent_event = 0
        cent_event_list = []
        high_mul, low_mul = get_centlimit(fpath, cent)
        #print (cent,high_mul,low_mul)
        multlist = np.loadtxt('./events/%s/mult_seta_max_%s.dat'%(comment,cpuid))
        for id, mult_tem in enumerate(multlist): 
            if mult_tem < high_mul and mult_tem > low_mul:
                cent_event_list.append(id)
                cent_event = cent_event+1

                print (cent_event,id,mult_tem)
            if cent_event >  nchoose:
                break

        cent_event_file.create_dataset('cent/{cent}/'.format(cent=cent),data=cent_event_list,dtype="i")
        eventlist=np.concatenate((eventlist, cent_event_list))
        print (centlist,cent,eventlist,cent_event_list)   
    cent_event_file.close()
    
    datalist=np.unique(eventlist)

    for eventid in datalist:
        source_path="event_{0:d}".format(int(eventid))
        target_path=source_path
        call("h5copy -i {0} -o {1} -s {2} -d {3}".format(all_events_path,cent_events_path,source_path,target_path),shell=True)
        print("h5copy -i {0} -o {1} -s {2} -d {3}".format(all_events_path,cent_events_path,source_path,target_path))


  
             
    




if __name__ == '__main__':
    start_time = time.time() 

    #first, generate many multiplicity data
    projectile,target = 'p','Pb'
    num_of_event=1

    #the following parameter all need to be confirmed further************************************************
    reduced_thickness=0
    beam_energy=8160
    cross_section=7.25

    bmin=0.0
    bmax=12.0
    grid_max= 10.0
    grid_step=0.10

    # grid_max= 9.9
    # grid_step=0.3
    
    eta_max=16
    eta_step=0.16
    
    mean_coeff=0.0
    std_coeff=2.9
    skew_coeff=6.0
    skew_type=1
    jacobian=0.75
    fluctuation=2.0
    
    nucleon_width=0.59
    constituent_width=0.3
    constituent_num=3
    min_dist=0.7
    
    cwd=os.getcwd()
    comment=sys.argv[1]
    if not os.path.exists(os.path.join('events',comment)): os.makedirs(os.path.join('events',comment))
    dir_num=10
    cpuid=0

    #if True:
    name=['projectile','target','num_of_event','reduced_thickness','beam_energy','cross_section',\
         'bmin','bmax','grid_max=','grid_step','eta_max','eta_step',\
         'mean_coeff','std_coeff','skew_coeff','skew_type','jacobian','fluctuation',\
         'nucleon_width','constituent_width','constituent_num','min_dist']
    value=[projectile,target,num_of_event,reduced_thickness,beam_energy,cross_section,\
         bmin,bmax,grid_max,grid_step,eta_max,eta_step,\
         mean_coeff,std_coeff,skew_coeff,skew_type,jacobian,fluctuation,\
         nucleon_width,constituent_width,constituent_num,min_dist]
        
    df=pd.Series(value,index=name)
    df.to_csv(os.path.join('events',comment,'parameters.csv'),sep=':')
                
    pool=mp.Pool(processes=dir_num)
    for threadid in range(dir_num):
        pool.apply_async(run_jet_position,args=(comment,cpuid,threadid,projectile,target,num_of_event,reduced_thickness, beam_energy,cross_section,\
                                        bmin,bmax,grid_max,grid_step,eta_max,eta_step,\
                                        mean_coeff,std_coeff,skew_coeff,skew_type,jacobian,fluctuation,\
                                        nucleon_width,constituent_width,constituent_num,min_dist))
    
    

    pool.close()
    pool.join()

    centlist=["0_1","1_5","0_5","5_10","10_20","20_40","40_60","60_80","80_100","0_100","20_50"]
    get_multi_list_and_centrality_cut2('.',dir_num,grid_max,grid_step,eta_max,eta_step,comment=comment,cpuid=cpuid)
    collect_event(centlist,comment=comment,cpuid=cpuid,nchoose=200)
    end_time = time.time()  
    print(f"time: {end_time - start_time} s")
