##def overlap2structure(
##    l_evecs1, l_evecs2,d_coordinates,l_rem,d_hessian,l_coordinates,
##    prefix,suffix,
##    ):
##
##    for mode in range(6,12):
##        evec1 = l_evecs1[mode]
##        evec2 = l_evecs2[mode]
##
##        lines = []
##
##        j = -1
##        for i in d_hessian.keys():
##            j += 1
##
##            if i in l_rem:
##                j -= 1
##
##            record = d_hessian[i]['record']
##            atom_no = d_hessian[i]['atom_no']
##            atom_name = d_hessian[i]['atom_name']
##            altloc = d_hessian[i]['altloc']
##            res_name = d_hessian[i]['res_name']
##            chain = d_hessian[i]['chain']
##            res_no = d_hessian[i]['res_no']
##            iCode = d_hessian[i]['iCode']
##            coordinate = d_hessian[i]['coordinate']
##            x = coordinate[0]
##            y = coordinate[1]
##            z = coordinate[2]
##            occupancy = 1.0
##            if i in l_rem:
##                bfactor = 99.99
##            else:
##                v1 = evec1[3*i:3*(i+1)]
##                v2 = evec2[3*j:3*(j+1)]
##                overlap = cosangle(v1,v2,)
##                bfactor = 50+50*overlap
##            element = d_hessian[i]['element']
##
##            line = '%6s%5i  %3s%s%3s %1s%4i%1s   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s  \n' %(
##                record.ljust(6),atom_no,
##                atom_name.ljust(3), altloc, res_name, chain, res_no, iCode,
##                x,y,z, occupancy, bfactor,element.rjust(2)
##                )
##
##            lines += [line]
##
##        fd = open('%s_resoverlap_mode%02i_%s.pdb' %(prefix,mode,suffix,), 'w')
##        fd.writelines(lines)
##        fd.close()
##
##    return    


def vmdarrow(eigenvectors, d_hessian, prefix, suffix,l_rem=[],):

    import math

    lines = ['draw color white\n']

    j = -1
    for i in d_hessian.keys():
        j += 1

        if i in l_rem:
            j -= 1
            continue

        coordinate = d_hessian[i]['coordinate']
        x = coordinate[0]
        y = coordinate[1]
        z = coordinate[2]
        vx = 20*eigenvectors[6][3*j+0]
        vy = 20*eigenvectors[6][3*j+1]
        vz = 20*eigenvectors[6][3*j+2]
        v_len = math.sqrt(vx**2+vy**2+vz**2)

        ## cylinder *and* cone
        if v_len > 1:
            v_cone = [
                20*eigenvectors[6][3*j+0]/v_len,
                20*eigenvectors[6][3*j+1]/v_len,
                20*eigenvectors[6][3*j+2]/v_len,
                ]
            v_cylinder = [
                20*eigenvectors[6][3*j+0]-v_cone[0],
                20*eigenvectors[6][3*j+1]-v_cone[1],
                20*eigenvectors[6][3*j+2]-v_cone[2],
                ]
        ## cone only
        else:
            v_cone = [
                20*eigenvectors[6][3*j+0],
                20*eigenvectors[6][3*j+1],
                20*eigenvectors[6][3*j+2],
                ]
            v_cylinder = [0,0,0,]

        ## cylinder
        if v_len > 1:
            line = 'draw cylinder {%f %f %f} {%f %f %f} radius 0.1\n' %(
                x,y,z,
                x+v_cylinder[0],y+v_cylinder[1],z+v_cylinder[2],
                )
            lines += [line]

        ## cone
        line = 'draw cone {%f %f %f} {%f %f %f} radius 0.15\n' %(
            x+v_cylinder[0],y+v_cylinder[1],z+v_cylinder[2],
            x+v_cylinder[0]+v_cone[0],y+v_cylinder[1]+v_cone[1],z+v_cylinder[2]+v_cone[2],
            )
        lines += [line]

    fd = open('%s_%s.vmd' %(prefix,suffix,), 'w')
    fd.writelines(lines)
    fd.close()

    return


def morph(eigenvectors, nframes, chains, d_coordinates, job, d_MODRES, l_rem = [], l_add = [], write_frames = False):

    if l_add != []:
        stop

    '''visualize the two extreme projections along a trajectory and interpolate n frames between them'''

    import math

    d_colors = {6:'red',7:'green',8:'blue',9:'yellow',10:'violet',11:'orange',}

    ## adjust for difference in length between eigenvectors
    dic = {}
    j = 0
    for i in range(len(eigenvectors)+3*len(l_rem)):
        if i in l_rem:
            dic[i] = 'Continue'
            j -= 1
        else:
            dic[i] = i+j


    for mode in range(6,12):

        eigenvector = eigenvectors[mode]

        ## pdb header
        output_vmd = []
        vmd_arrows = ['draw color %s\n' %(d_colors[mode])]

        ## loop over frames
        for frame in range(nframes):
            
##            output_frame = []
            
            ## frame header
##                output_vmd.append('HEADER    frame t= %4.3f\nMODEL        0\n' %(frame))
            output_vmd.append('HEADER    frame t= %4.3f\nMODEL     %4i\n' %(frame,frame+1))
            ## loop over coordinates
            i = 0
            chains = d_coordinates['chains'].keys()
            chains.sort()
            for chain in chains:
                ## assume sequential numbering of residues
                res_nos = d_coordinates['chains'][chain]['residues'].keys()
                res_nos.sort()
                for res_no in res_nos:
                    ## assume sequential iCodes
                    iCodes = d_coordinates['chains'][chain]['residues'][res_no]['iCodes'].keys()
                    iCodes.sort()
                    for iCode in iCodes:
                        altloc_residue = min(d_coordinates['chains'][chain]['residues'][res_no]['iCodes'][iCode]['altlocs'].keys())
                        res_name = d_coordinates['chains'][chain]['residues'][res_no]['iCodes'][iCode]['altlocs'][altloc_residue]['res_name']
                        atom_names = d_coordinates['chains'][chain]['residues'][res_no]['iCodes'][iCode]['res_names'][res_name]['atoms'].keys()
                        atom_names.sort()
                        for atom_name in atom_names:
                            altloc = min(d_coordinates['chains'][chain]['residues'][res_no]['iCodes'][iCode]['res_names'][res_name]['atoms'][atom_name]['altlocs'].keys())
                            coordinate = d_coordinates['chains'][chain]['residues'][res_no]['iCodes'][iCode]['res_names'][res_name]['atoms'][atom_name]['altlocs'][altloc]['coordinate']
                            atom_no = d_coordinates['chains'][chain]['residues'][res_no]['iCodes'][iCode]['res_names'][res_name]['atoms'][atom_name]['altlocs'][altloc]['atom_no']

                            x1 = coordinate[0]
                            y1 = coordinate[1]
                            z1 = coordinate[2]

                            ## eigenvector?
                            if 'i' in d_coordinates['chains'][chain]['residues'][res_no]['iCodes'][iCode]['res_names'][res_name]['atoms'][atom_name]['altlocs'][altloc].keys():
                                i = d_coordinates['chains'][chain]['residues'][res_no]['iCodes'][iCode]['res_names'][res_name]['atoms'][atom_name]['altlocs'][altloc]['i']
                                j = dic[i]
                                if j == 'Continue':
                                    continue
                                else:
                                    i=j
                                if i in l_add:
                                    stop
                                vx = 32*eigenvector[3*i+0]
                                vy = 32*eigenvector[3*i+1]
                                vz = 32*eigenvector[3*i+2]
                                sqlength = vx**2+vy**2+vz**2
                                len_v = math.sqrt(sqlength)
                                x2 = x1+(1-2*float(frame)/float(nframes))*vx
                                y2 = y1+(1-2*float(frame)/float(nframes))*vy
                                z2 = z1+(1-2*float(frame)/float(nframes))*vz
##                                    x += (float(frame)/float(nframes))*vx
##                                    y += (float(frame)/float(nframes))*vy
##                                    z += (float(frame)/float(nframes))*vz
                            else:
                                continue

                            ## MODRES?
                            MODRES = False
                            if chain in d_MODRES.keys():
                                if res_no in d_MODRES[chain].keys():
                                    if iCode in d_MODRES[chain][res_no].keys():
                                        if res_name != d_MODRES[chain][res_no][iCode]:
                                            print res_name, d_MODRES[chain][res_no][iCode]
                                            notexpected
                                        MODRES = True

                            ## ATOM or HETATM ?
                            if MODRES == False and res_name not in ['GLY','ALA','VAL','LEU','ILE','SER','THR','CYS','MET','ASP','GLU','ASN','GLN','HIS','LYS','ARG','PRO','PHE','TRP','TYR',]:
                                record = 'HETATM'
                            elif MODRES == True:
                                record = 'ATOM'
## add a flag to know if std res is atom or hetatm
                            elif MODRES == False and res_name in ['GLY','ALA','VAL','LEU','ILE','SER','THR','CYS','MET','ASP','GLU','ASN','GLN','HIS','LYS','ARG','PRO','PHE','TRP','TYR',]:
                                record = 'ATOM'

                            ## append atom line
                            line = '%6s%5i  %3s%s%3s %1s%4i%1s   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s  \n' %(
                                record.ljust(6),atom_no,
                                atom_name.ljust(3), altloc, res_name, chain, res_no, iCode,
                                x2,y2,z2, 1.0, sqlength,atom_name[0].rjust(2)
                                )
                            output_vmd += [line]

##                            for atom_name in ['N','C','CA',]:
##                                altloc = min(d_coordinates['chains'][chain]['residues'][res_no]['iCodes'][iCode]['res_names'][res_name]['atoms'][atom_name]['altlocs'].keys())
##                                atom_no = d_coordinates['chains'][chain]['residues'][res_no]['iCodes'][iCode]['res_names'][res_name]['atoms'][atom_name]['altlocs'][altloc]['atom_no']
##                                coordinate = d_coordinates['chains'][chain]['residues'][res_no]['iCodes'][iCode]['res_names'][res_name]['atoms'][atom_name]['altlocs'][altloc]['coordinate']
##                                x1 = coordinate[0]
##                                y1 = coordinate[1]
##                                z1 = coordinate[2]
##                                vx = 20*eigenvector[3*i+0]
##                                vy = 20*eigenvector[3*i+1]
##                                vz = 20*eigenvector[3*i+2]
##                                sqlength = vx**2+vy**2+vz**2
##                                len_v = math.sqrt(sqlength)
##                                x2 = x1+(1-2*float(frame)/float(nframes))*vx
##                                y2 = y1+(1-2*float(frame)/float(nframes))*vy
##                                z2 = z1+(1-2*float(frame)/float(nframes))*vz
##                                line = '%6s%5i  %3s%s%3s %1s%4i%1s   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s  \n' %(
##                                    record.ljust(6),atom_no,
##                                    atom_name.ljust(3), altloc, res_name, chain, res_no, iCode,
##                                    x2,y2,z2, 1.0, sqlength,atom_name[0].rjust(2)
##                                    )
##                                output_frame += [line]

                            if frame == 0:
                                if sqlength < 1:
                                    v_cone = [
                                        20*eigenvector[3*i+0],
                                        20*eigenvector[3*i+1],
                                        20*eigenvector[3*i+2],
                                        ]
                                    v_cylinder = [0,0,0,]
                                else:
                                    v_cone = [
                                        20*eigenvector[3*i+0]/len_v,
                                        20*eigenvector[3*i+1]/len_v,
                                        20*eigenvector[3*i+2]/len_v,
                                        ]
                                    v_cylinder = [
                                        20*eigenvector[3*i+0]-v_cone[0],
                                        20*eigenvector[3*i+1]-v_cone[1],
                                        20*eigenvector[3*i+2]-v_cone[2],
                                        ]
                                if sqlength > 1:
                                    line = 'draw cylinder {%f %f %f} {%f %f %f} radius 0.1\n' %(
                                        x1,y1,z1,
                                        x1+v_cylinder[0],y1+v_cylinder[1],z1+v_cylinder[2],
                                        )
                                    vmd_arrows += [line]
                                line = 'draw cone {%f %f %f} {%f %f %f} radius 0.15\n' %(
                                    x1+v_cylinder[0],y1+v_cylinder[1],z1+v_cylinder[2],
                                    x1+v_cylinder[0]+v_cone[0],y1+v_cylinder[1]+v_cone[1],z1+v_cylinder[2]+v_cone[2],
                                    )
                                vmd_arrows += [line]

                            i += 1

            output_vmd.append('TER\nENDMDL\n')

##            if mode == 7:
##                fd = open('%s.pdb' %(frame),'w')
##                fd.writelines(output_frame)
##                fd.close()

        fd = open(job+'_mode'+str(mode+1).zfill(2)+'.pdb', 'w')
        fd.writelines(output_vmd)
        fd.close()
        fd = open('%s_vmdarrow_mode%02i.src' %(job, mode+1), 'w')
        fd.writelines(vmd_arrows)
        fd.close()

    return


def gnuplot_histogram(d_data, fileprefix):

    import math, os

    l_xtics = d_data.keys()
    l_xtics.sort()

    ##
    ## write data
    ##
    yrange = []
    gnuplotdata = []
    for i in range(len(l_xtics)):
        res_name = l_xtics[i]
        for y in d_data[res_name]:
            gnuplotdata += ['%f %f\n' %(float(i),y)]
            yrange += [y]
    ymin = min(yrange)
    ymax = max(yrange)
    fd = open('gnuplot.data','w')
    fd.writelines(gnuplotdata)
    fd.close()

    ##
    ## calculate statistics
    ##
    gnuplot_statistics = []
    for i in range(len(l_xtics)):
        res_name = l_xtics[i]
        n = len(d_data[res_name])
        if n <= 1:
            continue
        sumx = 0
        sumxx = 0
        for x in d_data[res_name]:
            sumx += x
            sumxx += x**2
        average = sumx/n
        SS = sumxx-(sumx**2)/n
        MSE = SS / (n-1)
        if MSE < 0:
            SE = 0 ## temp!!! check the equation!!!
        else:
            SE = math.sqrt(MSE/n)
        gnuplot_statistics += ['%f %f %f\n' %(float(i),average,SE)]
    ## write statistics
    fd = open('gnuplot.statistics','w')
    fd.writelines(gnuplot_statistics)
    fd.close()

    ##
    ## write gnuplot settings
    ##
    gnuplotsettings = []
    gnuplotsettings += [
        'set terminal postscript eps enhanced color "Helvetica" 48\n',
        'set output "gnuplot.ps"\n',
        'set size 4,4\n', ## scale 400%
        'set autoscale fix\n', ## scale axes to include min and max *only* and *not* the next tic
        'set xlabel "classification"\n',
        'set ylabel "max overlap"\n',
        ]
    line_xtic = 'set xtics ('
    for xtic in l_xtics:
        line_xtic += '"%s" %s, ' %(xtic, l_xtics.index(xtic))
    line_xtic = line_xtic[:-2]+')\n'
    gnuplotsettings += [
        line_xtic,
        'set xtics rotate\n',
    ]
    gnuplotsettings += [
        'plot [-1:%i][%f:%f] "gnuplot.data" lt 0 ps 2 pt 2' %(len(l_xtics)+1, ymin, ymax),
        ', "gnuplot.statistics" lt 1 lc 0 ps 0 pt 0 w errorb\n',
        ]
    ## write gnuplot settings
    fd = open('gnuplot.settings','w')
    fd.writelines(gnuplotsettings)
    fd.close()

    ##
    ## execute gnuplot settings
    ##
    os.system('/software/bin/gnuplot gnuplot.settings')
    ## convert postscript to portable network graphics
    os.system('convert gnuplot.ps %s.png' %(fileprefix))
    os.remove('gnuplot.ps')
    os.remove('gnuplot.data')
    os.remove('gnuplot.settings')
    os.remove('gnuplot.statistics')

    return


def CA2CB(remCA,d_coordinates,d_hessian):
    chain = d_hessian[remCA]['chain']
    res_no = d_hessian[remCA]['res_no']
    iCode = d_hessian[remCA]['iCode']
    res_name = d_hessian[remCA]['res_name']
    l_rem = []
    for atom_name in d_coordinates['chains'][chain]['residues'][res_no]['iCodes'][iCode]['res_names'][res_name]['atoms'].keys():
        altloc = min(d_coordinates['chains'][chain]['residues'][res_no]['iCodes'][iCode]['res_names'][res_name]['atoms'][atom_name]['altlocs'].keys())
        if 'i' in d_coordinates['chains'][chain]['residues'][res_no]['iCodes'][iCode]['res_names'][res_name]['atoms'][atom_name]['altlocs'][altloc].keys():
            i = d_coordinates['chains'][chain]['residues'][res_no]['iCodes'][iCode]['res_names'][res_name]['atoms'][atom_name]['altlocs'][altloc]['i']
            l_rem += [i]
    return l_rem


def parse_pdb(lines, parse_chains, parse_remarks = True):

    import Numeric, sets

    d_coordinates = {'chains':{}} ## ATOM, HETATM, MODEL
    d_REMARK350 = {} 
    d_secondary = {'HELIX':{},'SHEET':{},} ## HELIX, SHEET
    d_primary = {'SEQRES':sets.Set(),'MODRES':{}}
    d_ligands = {'chains':{}}
    d_atomnos = {}
    accept_missing_residues = False
    accept_altloc = False

    for i in range(len(lines)):

        line = lines[i]

        record = line[:6].strip()

        if record == 'ATOM':
            altloc = line[16]
            if altloc != ' ' and accept_altloc == True:
                print 'Alternative locations of atoms exist. Do you want to proceed (Yes/No)?'
                proceed = raw_input()
                if proceed.upper() in ['YES','Y']:
                    accept_altloc = True
                elif proceed.upper() in ['NO','N']:
                    accept_altloc = False
                    raise
                else:
                    raise 'invalid answer'

            d_coordinates = parse_ATOM(line, d_coordinates)
            continue

        elif record == 'REMARK' and parse_remarks == True:

            remark = int(line[6:10].strip())

            if remark == 350:

                d_REMARK350 = parse_REMARK350(lines,i, d_REMARK350)

            elif remark == 465:

                if accept_missing_residues == False:
                    print 'Residues are missing (i.e. REMARK465 records are present). Do you want to proceed (Yes/No)?'
                    proceed = raw_input()
                    if proceed.upper() in ['YES','Y']:
                        accept_missing_residues = True
                    elif proceed.upper() in ['NO','N']:
                        accept_missing_residues = False
                        raise
                    else:
                        raise 'invalid answer'

        elif record == 'SEQRES':
            chain = line[11]
            d_primary['SEQRES'] |= sets.Set(chain)
                                
        elif record == 'MODRES':
            chain = line[16]
            res_no = int(line[18:22])
            iCode = line[22]
            res_name = line[12:15]
            if chain not in d_primary['MODRES'].keys():
                d_primary['MODRES'][chain] = {}
            if res_no not in d_primary['MODRES'][chain].keys():
                d_primary['MODRES'][chain][res_no] = {}
            if iCode not in d_primary['MODRES'][chain][res_no].keys():
                d_primary['MODRES'][chain][res_no][iCode] = res_name

        elif record == 'HELIX':
            chain = line[19]
            if not d_secondary['HELIX'].has_key(chain):
                d_secondary['HELIX'][chain] = {'range':[],'list':[]}
            d_secondary['HELIX'][chain]['range'] += [[int(line[21:25]), int(line[33:37])]]
            d_secondary['HELIX'][chain]['list'] += range(int(line[21:25]), int(line[33:37])+1)
            continue

        elif record == 'SHEET':
            chain = line[21]
            if not d_secondary['SHEET'].has_key(chain):
                d_secondary['SHEET'][chain] = {'range':[]} 
            d_secondary['SHEET'][chain]['range'] += [[int(line[22:26]), int(line[33:37])]]
            continue

        elif record == 'MODEL':
            model = int(line[10:14])
            continue

        elif record == 'HETATM':
            MODRES = False
            chain = line[21]
            res_no = int(line[22:26].strip())
            iCode = line[26]
            res_name = line[17:20].strip()
            if res_name in ['HOH','DOD',]:
                continue
            if chain in d_primary['MODRES'].keys():
                if res_no in d_primary['MODRES'][chain].keys():
                    if iCode in d_primary['MODRES'][chain][res_no].keys():
                        if res_name != d_primary['MODRES'][chain][res_no][iCode]:
                            print res_name, d_primary['MODRES'][chain][res_no][iCode]
                            notexpected
                        MODRES = True
            if MODRES == True:
                d_coordinates = parse_ATOM(line, d_coordinates)
            else:
                d_ligands = parse_ATOM(line,d_ligands)

            continue

    return d_REMARK350, d_primary, d_secondary, d_coordinates, d_ligands


def parse_REMARK350(lines, i, d_REMARK350):

    import sets

    line = lines[i]

    if line[11:23] == 'BIOMOLECULE:':
        biomolecules = line[23:80].replace(' ','').split(',')

        chains = sets.Set()

        for j in range(i+1,len(lines)):

            if lines[j][:10] != 'REMARK 350' or lines[j][11:23] == 'BIOMOLECULE:':
                break

            elif lines[j][11:41] == 'APPLY THE FOLLOWING TO CHAINS:':
                chains = sets.Set() ## e.g. 1upp.pdb (vs 8ruc.pdb)
                line_chains = lines[j][41:80]
                chains |= parse_REMARK350_chains(line_chains)

            elif lines[j][11:53] == 'IN ADDITION APPLY THE FOLLOWING TO CHAINS:':
                line_chains = lines[j][53:80]
                chains |= parse_REMARK350_chains(line_chains)

            elif ',' in lines[j][11:80]:
                if 'APPLY THE FOLLOWING TO CHAINS:' in [lines[j-1][11:41],lines[j-2][11:41]]:
                    line_chains = lines[j][11:80]
                    chains |= parse_REMARK350_chains(line_chains)

            ## count and parse chain transformations
            ## accept SMTRY3 to accomodate for pdb entries from 1996 and prior (e.g. 1thj.pdb)
            elif lines[j][13:19] in ['BIOMT3','SMTRY3']:

                matrixno = int(lines[j][19:24])
                ## parse transformation matrix
                matrixrow1 = lines[j-2][24:].split()
                matrixrow2 = lines[j-1][24:].split()
                matrixrow3 = lines[j-0][24:].split()
                matrixrows = [matrixrow1,matrixrow2,matrixrow3,]
##                ## find out whether transformation matrix yields a transformation
##                transformation = False
##                for k in range(3):
##                    ## add a zero translation vector if a translation vector is not given
##                    if len(matrixrows[k]) == 3:
##                        matrixrows[k] += [0.]
##                    if float(matrixrows[k][k]) == 1. and float(matrixrows[k][3]) == 0.:
##                        continue
##                    else:
##                        transformation = True

                ## append transformation matrix to dictionary
                for biomolecule in biomolecules:

                    biomolecule = int(biomolecule)

                    ## biomolecule
                    if biomolecule not in d_REMARK350.keys():
                        d_REMARK350[biomolecule] = {}

                    ## biomolecule > matrices
                    if 'matrices' not in d_REMARK350[biomolecule].keys():
                        d_REMARK350[biomolecule]['matrices'] = {}
                    ## matrices > matrixno > matrix
                    d_REMARK350[biomolecule]['matrices'][matrixno] = matrixrows

                    ## biomolecule > chains
                    if 'chains' not in d_REMARK350[biomolecule].keys():
                        d_REMARK350[biomolecule]['chains'] = {}
                    for chain in chains:
                        ## chains > chain
                        if chain not in d_REMARK350[biomolecule]['chains'].keys():
                            d_REMARK350[biomolecule]['chains'][chain] = sets.Set()
                        d_REMARK350[biomolecule]['chains'][chain] |= sets.Set([matrixno])


    return d_REMARK350


def parse_REMARK350_chains(line_chains):

    import sets

    ## if sentence necessary due to e.g. 1qgc
    ## REMARK 350 APPLY THE FOLLOWING TO CHAINS: 1 2 3 4 5
    if ',' not in line_chains:
        chains = line_chains.split()
    else:
        ## remove 'AND' from the line of chains (e.g. problem with 1rhi)
        ## replace '.' in the line of chains (e.g. problem with 1rbo and 1qgc)
        chains = line_chains.replace('AND',',').replace('.',',').replace(' ',',').split(',')

    ## loop removal of blank chains necessary due to e.g. 2g8g
    ## REMARK 350 APPLY THE FOLLOWING TO CHAINS: A, ,
    for x in range(100):
        if '' in chains:
            chains.remove('')
        else:
            break

    for j in range(len(chains)):
        chain = chains[j]
        if chain == 'NULL':
            chains[j] = ' '

    return sets.Set(chains)


def loop_and_identify_biomolecules(i, lines):

    for j in range(i-1,-1,-1):

        if lines[j][:10] != 'REMARK 350':
            return False
        if lines[j][11:23] == 'BIOMOLECULE:':
            return True


def parse_ATOM(line, d_coordinates):

    import Numeric

    record = line[:6].strip()
    atom_no = int(line[6:11])
    atom_name = line[12:16].strip()
    altloc = line[16]
    res_name = line[17:20].strip()
    chain = line[21]
    res_no = int(line[22:26])
    iCode = line[26]
    x = float(line[30:38])
    y = float(line[38:46])
    z = float(line[46:54])
    occupancy = float(line[54:60])
    beta = float(line[60:66])
    if altloc != ' ':
        print 'alternate atom locations for atom %s of residue %s%s%s in chain %s' %(atom_name,res_name,res_no,iCode,chain,)
    if altloc not in [' ','A',]:
        return d_coordinates
    
    element = line[76:78].strip()

    coordinate = Numeric.array([x, y, z])

    if not chain in d_coordinates['chains'].keys():
        d_coordinates['chains'][chain] = {'residues':{}}
    if not res_no in d_coordinates['chains'][chain]['residues'].keys():
        d_coordinates['chains'][chain]['residues'][res_no] = {'iCodes':{}}
    if not iCode in d_coordinates['chains'][chain]['residues'][res_no]['iCodes'].keys():
        d_coordinates['chains'][chain]['residues'][res_no]['iCodes'][iCode] = {'res_names':{},'altlocs':{}}
    if not res_name in d_coordinates['chains'][chain]['residues'][res_no]['iCodes'][iCode]['res_names'].keys():
        d_coordinates['chains'][chain]['residues'][res_no]['iCodes'][iCode]['res_names'][res_name] = {'atoms':{}}
    if not atom_name in d_coordinates['chains'][chain]['residues'][res_no]['iCodes'][iCode]['res_names'][res_name]['atoms'].keys():
        d_coordinates['chains'][chain]['residues'][res_no]['iCodes'][iCode]['res_names'][res_name]['atoms'][atom_name] = {'altlocs':{}}

    d_coordinates['chains'][chain]['residues'][res_no]['iCodes'][iCode]['altlocs'][altloc] = {
        'res_name':res_name
        }
    d_coordinates['chains'][chain]['residues'][res_no]['iCodes'][iCode]['res_names'][res_name]['atoms'][atom_name]['altlocs'][altloc] = {
        'coordinate':coordinate,'element':element,'atom_no':atom_no,
        'occupancy':occupancy,'beta':beta,'record':record,
        }

    return d_coordinates


def align_vectors_of_different_length(eigenvectors_nonperturbed, eigenvectors_perturbed, l_rem = [], l_add = []):

    stop

    if len(eigenvectors_nonperturbed)+3*len(l_add) != len(eigenvectors_perturbed)+3*len(l_rem):
        print len(eigenvectors_nonperturbed), len(eigenvectors_perturbed), 3*len(l_rem)
        stop

    l_eigenvectors_nonperturbed_aligned = []
    for mode_nonperturbed in range(len(eigenvectors_nonperturbed)):
        vector_nonperturbed = list(eigenvectors_nonperturbed[mode_nonperturbed])
        l_rem.sort()
        l_rem.reverse()
        for i in l_rem:
            for j in range(3-1,-1,-1):
                del vector_nonperturbed[3*i+j]
        l_eigenvectors_nonperturbed_aligned += [vector_nonperturbed]

    l_eigenvectors_perturbed_aligned = []
    for mode_perturbed in range(len(eigenvectors_perturbed)):
        vector_perturbed = list(eigenvectors_perturbed[mode_perturbed])
        l_add.sort()
        l_add.reverse()
        for i in l_add:
            for j in range(3-1,-1,-1):
                del vector_perturbed[3*i+j]
        l_eigenvectors_perturbed_aligned += [vector_perturbed]

    return l_eigenvectors_nonperturbed_aligned, l_eigenvectors_perturbed_aligned


def overlap_calculation(eigenvectors_perturbed, eigenvectors_nonperturbed, eigenvalues_perturbed, eigenvalues_nonperturbed, l_rem = [], l_add = []):

    ## calculate fewer overlaps to speed up things. all overlaps only calculated to determine mode of max overlap

    import math, Numeric

    ##
    ## reset lists
    ##
    overlaps = []
    max_overlaps = []
    perturbed_modes_of_max_overlap = []
    delta_perturbed_eigenvalues_of_max_overlap = []

    ##
    ## combinatorial loop over modes
    ##
    moderanges_nonperturbed = range(6,12)+[len(eigenvectors_nonperturbed)-3*len(l_rem)-1]
    for mode_nonperturbed in range(len(eigenvectors_nonperturbed)):
        max_overlap = 0
        overlaps_per_mode_perturbed = []
        vector_nonperturbed = list(eigenvectors_nonperturbed[mode_nonperturbed])
        ## adjust for different lengths of eigenvectors
        l_rem.sort()
        l_rem.reverse()
        for i in l_rem:
            for j in range(3-1,-1,-1):
                del vector_nonperturbed[3*i+j]
        overlaps_per_mode_perturbed = [0 for mode_perturbed in range(len(eigenvectors_perturbed))]
        moderanges_perturbed = range(6,12)+[len(eigenvectors_nonperturbed)-3*len(l_rem)-1]
        for mode_perturbed in range(len(eigenvectors_perturbed)):
            vector_perturbed = list(eigenvectors_perturbed[mode_perturbed])
            ## adjust for different lengths of eigenvectors
            l_add.sort()
            l_add.reverse()
            for i in l_add:
                for j in range(3-1,-1,-1):
                    del vector_perturbed[3*i+j]
            if mode_nonperturbed in moderanges_nonperturbed and mode_perturbed in moderanges_perturbed:
                overlap = math.fabs(cosangle(vector_nonperturbed, vector_perturbed))
            ## reduce moderanges to speed up things significantly
            else:
                overlap = 0
            overlaps_per_mode_perturbed[mode_perturbed] = overlap
            if mode_nonperturbed == mode_perturbed:
                overlaps.append(overlap)
        ## identify max overlap in list of overlaps
        max_overlap = max(overlaps_per_mode_perturbed)
        max_overlaps.append(max_overlap)
        ## identify mode of max overlap
        perturbed_mode_of_max_overlap = overlaps_per_mode_perturbed.index(max_overlap)
        perturbed_modes_of_max_overlap.append(perturbed_mode_of_max_overlap+1)
        ## identify eigenvalue of mode of max overlap
        delta_perturbed_eigenvalue_of_max_overlap = eigenvalues_perturbed[perturbed_mode_of_max_overlap]-eigenvalues_nonperturbed[mode_nonperturbed]
        delta_perturbed_eigenvalues_of_max_overlap.append(delta_perturbed_eigenvalue_of_max_overlap)

    return overlaps, max_overlaps, perturbed_modes_of_max_overlap, delta_perturbed_eigenvalues_of_max_overlap


def sigmoid(x, cutoff, slope=1):

    import math

    y = 1. / ( 1. + math.exp( slope*(x-cutoff) ) )

    return y


def platonicsolid(origo, polyhedron = None):

    import math

    a = 1 ## edge length

    if polyhedron == 'tetrahedron':
##        v1 = origo+[-1./2,-1./(2*math.sqrt(3)),-1./(2*math.sqrt(6))]
##        v2 = origo+[ 1./2,-1./(2*math.sqrt(3)),-1./(2*math.sqrt(6))]
##        v3 = origo+[ 0./2, 1./(1*math.sqrt(3)),-1./(2*math.sqrt(6))]
##        v4 = origo+[ 0./2, 0./(2*math.sqrt(3)), 3./(2*math.sqrt(6))]
        v1 = origo+[+1,+1,+1,]
        v2 = origo+[-1,-1,+1,]
        v3 = origo+[-1,+1,-1,]
        v4 = origo+[+1,-1,-1,]
        coordinates = [v1,v2,v3,v4,]
    elif polyhedron == 'hexahedron':
        c1 = origo+[+1*a,+1*a,+1*a,]
        c2 = origo+[+1*a,+1*a,-1*a,]
        c3 = origo+[+1*a,-1*a,+1*a,]
        c4 = origo+[+1*a,-1*a,-1*a,]
        c5 = origo+[-1*a,+1*a,+1*a,]
        c6 = origo+[-1*a,+1*a,-1*a,]
        c7 = origo+[-1*a,-1*a,+1*a,]
        c8 = origo+[-1*a,-1*a,-1*a,]
        coordinates = [c1,c2,c3,c4,c5,c6,c7,c8,]
    elif polyhedron == 'octahedron':
        c1 = origo+[+1*a,0,0,]
        c2 = origo+[-1*a,0,0,]
        c3 = origo+[0,+1*a,0,]
        c4 = origo+[0,-1*a,0,]
        c5 = origo+[0,0,+1*a,]
        c6 = origo+[0,0,-1*a,]
        coordinates = [c1,c2,c3,c4,c5,c6,]
    elif polyhedron == 'icosahedron':
        phi = (1+math.sqrt(5))/2.## golden ratio
        c1 = origo+[0,+1,+phi,]
        c2 = origo+[0,+1,-phi,]
        c3 = origo+[0,-1,+phi,]
        c4 = origo+[0,-1,-phi,]
        c5 = origo+[+1,+phi,0,]
        c6 = origo+[+1,-phi,0,]
        c7 = origo+[-1,+phi,0,]
        c8 = origo+[-1,-phi,0,]
        c9 = origo+[+phi,0,+1,]
        c10 = origo+[+phi,0,-1,]
        c11 = origo+[-phi,0,+1,]
        c12 = origo+[-phi,0,-1,]
        coordinates = [c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,]
    elif polyhedron == 'dodecahedron':
        vertices = 20
        stop_too_many_vertices

    return coordinates
    

def eigenv_calccomb(matrix_hessian, jobid, verbose, l_rem = []):

    '''Calculates eigenvectors and eigenvalues of a matrix.'''
    if verbose == True:
        print 'calculating eigenvalues and eigenvectors of the '+str(len(matrix_hessian))+'x'+str(len(matrix_hessian))+' Hessian matrix'
    import math

##    import time
##
##    import numpy
##    
##    t1 = time.clock()
##    eigen_tuple = numpy.linalg.eig(matrix_hessian)
##    t2 = time.clock()
##    print t2-t1
##    eigenvalues = list(eigen_tuple[0])
##    eigenvectors = list(eigen_tuple[1])
##    eigenvectors = numpy.transpose(eigenvectors) ## transpose when diagonalization with numpy instead of Numeric
##    eigen_list = zip(eigenvalues, eigenvectors)
##    eigen_list.sort()
##    eigenvalues = [eigen_list[eigen][0] for eigen in range(len(eigen_list))]
##    eigenvectors = [eigen_list[eigen][1] for eigen in range(len(eigen_list))]
##
##    t1 = time.clock()
##    eigen_tupleh = numpy.linalg.eigh(matrix_hessian)
##    t2 = time.clock()
##    print t2-t1
##    eigenvaluesh = list(eigen_tupleh[0])
##    eigenvectorsh = list(eigen_tupleh[1])
##    eigenvectors = numpy.transpose(eigenvectors) ## transpose when diagonalization with numpy instead of Numeric
##    eigen_listh = zip(eigenvaluesh, eigenvectorsh)
##    eigen_listh.sort()
##    eigenvaluesh = [eigen_listh[eigen][0] for eigen in range(len(eigen_listh))]
##    eigenvectorsh = [eigen_listh[eigen][1] for eigen in range(len(eigen_listh))]
##
####    import Numeric, LinearAlgebra
####
####    t1 = time.clock()
####    eigen_tuple = LinearAlgebra.eigenvectors(matrix_hessian)
####    t2 = time.clock()
####    print t2-t1
####    eigenvalues = list(eigen_tuple[0])
####    eigenvectors = list(eigen_tuple[1])
######    eigenvectors = numpy.transpose(eigenvectors) ## transpose when diagonalization with numpy instead of Numeric
####    eigen_list = zip(eigenvalues, eigenvectors)
####    eigen_list.sort()
####    eigenvalues = [eigen_list[eigen][0] for eigen in range(len(eigen_list))]
####    eigenvectors = [eigen_list[eigen][1] for eigen in range(len(eigen_list))]
####
####    t1 = time.clock()
####    eigen_tupleh = LinearAlgebra.Heigenvectors(matrix_hessian)
####    t2 = time.clock()
####    print t2-t1
####    eigenvaluesh = list(eigen_tupleh[0])
####    eigenvectorsh = list(eigen_tupleh[1])
######    eigenvectors = numpy.transpose(eigenvectors) ## transpose when diagonalization with numpy instead of Numeric
####    eigen_listh = zip(eigenvaluesh, eigenvectorsh)
####    eigen_listh.sort()
####    eigenvaluesh = [eigen_listh[eigen][0] for eigen in range(len(eigen_listh))]
####    eigenvectorsh = [eigen_listh[eigen][1] for eigen in range(len(eigen_listh))]
##
##    print eigenvalues[0], eigenvalues[-1]
##    print eigenvaluesh[0], eigenvaluesh[-1]
##
##    print eigenvectors[0][0], eigenvectors[0][-1], eigenvectors[-1][0], eigenvectors[-1][-1]
##    print eigenvectorsh[0][0], eigenvectorsh[0][-1], eigenvectorsh[-1][0], eigenvectorsh[-1][-1]
##
##    print eigen_tuple[0]
##    print min(eigen_tuple[0])
##    print max(eigen_tuple[0])
##
##    stop

    import numpy
    ## diagonalize hessian matrix
    eigen_tuple = numpy.linalg.eigh(matrix_hessian)
    ## parse eigenvalues and eigenvectors
    eigenvalues = list(eigen_tuple[0])
    eigenvectors = list(eigen_tuple[1])
    eigenvectors = numpy.transpose(eigenvectors) ## transpose when diagonalization with numpy instead of Numeric
    ## organize eigenvalues and eigenvectors in list
    eigen_list = zip(eigenvalues, eigenvectors)
    ## sort list
    eigen_list.sort()
    ## parse sorted eigenvalues and eigenvectors
    eigenvalues = [eigen_list[eigen][0] for eigen in range(len(eigen_list))]
    eigenvectors = [eigen_list[eigen][1] for eigen in range(len(eigen_list))]
    if verbose == True:
        lines = ['rows=modes, cols=coordinates\n']
        for mode in range(6,len(eigenvectors)):
            lines += [str(eigenvectors[mode])+'\n']
        fd = open('%s_eigenvectors.txt' %(jobid),'w')
        fd.writelines(lines)
        fd.close()

    ## calculate length of mode 7 (equals 1 when using module linearalgebra)
    len7 = math.sqrt(sum(numpy.array(eigenvectors[6])*numpy.array(eigenvectors[6])))

    ## loop over modes
    for i in range(7+3*len(l_rem),len(eigenvalues)-1):

        ## calculate length of mode i
        leni = math.sqrt(sum(numpy.array(eigenvectors[i])*numpy.array(eigenvectors[i])))

        ## scale length of mode i relative to length of mode 7
        try:
            lenfactor = (len7/leni)/(eigenvalues[i]/eigenvalues[6+3*len(l_rem)])
        except:
            print eigenvalues
            print l_rem
            print i, len7, leni, eigenvalues[i], eigenvalues[6+3*len(l_rem)]
            stop_float_divison
        for j in range(len(eigenvectors[i])):
            eigenvectors[i][j] *= lenfactor

    ## copy lists of eigenvectors to eigenvectors_combined
    eigenvectors_combined = []
    for mode in range(len(eigenvalues)):
        eigenvectors_combined.append(list(eigenvectors[mode]))
    ## change mode i to be the sum of modes 7 to i
    for mode in range(7,len(eigenvalues)):
        for coordinate in range(len(eigenvalues)):
            eigenvectors_combined[mode][coordinate] += eigenvectors_combined[mode-1][coordinate]

    if verbose == True:
        print 'calculated eigenvalues and eigenvectors of the '+str(len(matrix_hessian))+'x'+str(len(matrix_hessian))+' Hessian matrix'

    return eigenvectors, eigenvalues, eigenvectors_combined


def datadic_return(N,range1=None,range2=None,):

    if range1 == None:
        range1 = range(N)
    if range2 == None:
        range2 = range(N)

    import Numeric

    datadic = {
        'emo': {'data': None, 'name': 'change in eigenvalue of max overlap', 'diagonal': 0, 'zrange': ['','']},
        'mmo': {'data': None, 'name': 'perturbed modes of max overlap', 'diagonal': 'mode', 'zrange': [7,12]},
        'overlaps_single': {'data': None, 'name': 'overlaps of single modes', 'diagonal': 'max', 'zrange': ['','']},
        'overlaps_combined': {'data': None, 'name': 'overlaps of combined modes', 'diagonal': 'max', 'zrange': ['','']},
        'overlaps_max': {'data': None, 'name': 'max overlaps of single modes', 'diagonal': 'max', 'zrange': ['','']},
        'eigenvalues_perturbed': {'data': None, 'name': 'eigenvalues of perturbed modes', 'diagonal': 'max', 'zrange': ['','']},
        }
    for key in datadic.keys():
        l = []
        for mode in range(3*N):
            if mode not in range(6,18)+[3*N-1]:
                l += ['']
            else:
                matrix = Numeric.zeros((N,N),typecode='d')
                for row in range1:
                    for col in range2:
                        matrix[row][col] = 0
                l += [matrix]
        datadic[key]['data'] = l
    
    return datadic


def multiple_cpus_write(datadic,filename,remres1):

    lines = []
    for key in datadic:
        lines.append('data '+key+'\n')
        for mode in range(6,12):
            lines.append('mode '+str(mode+1)+'\n')
            for remres2 in range(len(datadic[key]['data'][mode][remres1])):
                lines.append(str(datadic[key]['data'][mode][remres1][remres2])+'\n')
        if key == 'overlaps_combined':
            mode = len(datadic[key]['data'])-1
            lines.append('mode '+str(mode+1)+'\n')
            for remres2 in range(len(datadic[key]['data'][mode][remres1])):
                lines.append(str(datadic[key]['data'][mode][remres1][remres2])+'\n')
    fd = open(filename, 'w')
    fd.writelines(lines)
    fd.close()

    return


def multiple_cpus_read(datadic,filename,remres1):

    import os

    if os.path.isfile(filename):
        ## read datafile to data variables before you continue
        fd = open(filename, 'r')
        lines = fd.readlines()
        fd.close()
        for i in range(len(lines)):
            line = lines[i]
            if line[:4] == 'data':
                key = line[5:-1]
            elif line[:4] == 'mode':
                line_number_data = i+1
                mode = int(line[5:-1])-1
            elif i >= line_number_data:
                remres2 = i - line_number_data
                if remres2 < remres1:
                    continue
                datadic[key]['data'][mode][remres1][remres2] = float(line)
                datadic[key]['data'][mode][remres2][remres1] = float(line)
            else:
                print line, i, line_number_data, remres1, remres2
                stop

        Continue = True

    ## create data files (necessary to create data files if multithreading/parallel processing dependent on existence of files to determine if calculation on certain residues has been performed)
    else:
        Continue = False
        fd = open(filename, 'w')
        fd.close()

    return datadic, Continue



def pdb_import(pdb, path_pdb):

    
    import os, urllib2

    ## local dir
    if os.path.isfile('%s.pdb' %(pdb)):
        print 'importing pdb %s from local dir' %pdb
        fd = open(pdb+'.pdb', 'r')
        lines = fd.readlines()
        fd.close()
    elif os.path.isfile('%s.pdb' %(pdb.lower())):
        print 'importing pdb %s from local dir' %pdb
        fd = open(pdb.lower()+'.pdb', 'r')
        lines = fd.readlines()
        fd.close()
    elif os.path.isfile('%s.pdb' %(pdb.upper())):
        print 'importing pdb %s from local dir' %pdb
        fd = open(pdb.upper()+'.pdb', 'r')
        lines = fd.readlines()
        fd.close()
    elif os.path.isfile('%s%s/pdb%s.ent' %(path_pdb, pdb.lower()[1:3], pdb.lower())):
        print 'importing pdb %s from local pdb repository' %pdb
        fd = open('%s%s/pdb%s.ent' %(path_pdb, pdb.lower()[1:3], pdb.lower()), 'r')
        lines = fd.readlines()
        fd.close()
    else:
        print 'importing pdb %s from www.pdb.org' %pdb
        url = urllib2.urlopen('http://www.pdb.org/pdb/files/%s.pdb' %(pdb))
        lines = url.readlines()

    return lines


def cosangle(v1, v2):

    ## Numeric arrays are not used, because they are slow!

    import math, Numeric

    if len(v1) != len(v2):
        print len(v1), len(v2)
        stop

    numerator = 0
    for i in range(len(v1)):
        numerator += v1[i]*v2[i]

    denominator1 = 0
    denominator2 = 0
    for i in range(len(v1)):
        denominator1 += v1[i]*v1[i]
        denominator2 += v2[i]*v2[i]
    denominator = math.sqrt(denominator1*denominator2)

    if denominator == 0: ## temp!!!
        cosang = 0.
        stop
    else:
        cosang = numerator / denominator

    return cosang


def topology(
    fileprefix, chains, d_secondary, axlen
    ):

    import os

    bl = [390,833] ## bottom left coordinate
    tr = [1050,173] ## top right coordinate
    s = 14 ## space between plot and topology
    w = (tr[0]-bl[0])/axlen ## width of squares

    lines = 'convert %s.png' %(fileprefix)
    for chain in chains:
        for secelm in d_secondary.keys():
            if secelm == 'SHEET':
                lines += ' -fill "rgb(0,0,255)" '
            elif secelm == 'HELIX':
                lines += ' -fill "rgb(255,0,0)" '
            else:
                continue
            if chain in d_secondary[secelm]:
                for secelmrange in d_secondary[secelm][chain]['range']:
                    ## bottom
                    lines += ' -draw "line %s,%s %s,%s" ' %(int(bl[0]+w*(secelmrange[0]-1)), bl[1]+s, int(bl[0]+w*secelmrange[1]), bl[1]+s)
                    ## top
                    lines += ' -draw "line %s,%s %s,%s" ' %(int(bl[0]+w*(secelmrange[0]-1)), tr[1]-s, int(bl[0]+w*secelmrange[1]), tr[1]-s)
                    ## left
                    lines += ' -draw "line %s,%s %s,%s" ' %(bl[0]-s, int(bl[1]-w*(secelmrange[0]-1)), bl[0]-s, int(bl[1]-w*secelmrange[1]))
                    ## right
                    lines += ' -draw "line %s,%s %s,%s" ' %(tr[0]+s, int(bl[1]-w*(secelmrange[0]-1)), tr[0]+s, int(bl[1]-w*secelmrange[1]))
    if lines != 'convert %s.png' %(fileprefix):
        lines += '%s.png' %(fileprefix)
        os.system(lines)
            
    return


def post_perturbation_plot(
    datadic, jobid, chains, cutoff_distance,
    d_hessian, N, d_secondary, d_coordinates,
    winsize = 1,
    ):

    import math

    ##
    ## delta overlap vs residue number
    ##
    d_overlaps = {}
    for mode in range(6,12):

        l_overlaps = []
        for i in range(len(datadic['overlaps_single']['data'][mode])):
##            if datadic['overlaps_single']['data'][mode][i][i] != 0:
##                stop
            overlap = sum(datadic['overlaps_single']['data'][mode][i])/(len(datadic['overlaps_single']['data'][mode][i])-1)
            l_overlaps += [overlap]

        min_eigval = 'N/A'
        max_eigval = 0
        for i in range(len(datadic['eigenvalues_perturbed']['data'][mode])):
            for j in range(len(datadic['eigenvalues_perturbed']['data'][mode][i])):
                if i == j:
                    continue
                eigval = datadic['eigenvalues_perturbed']['data'][mode][i][j]
                if eigval < min_eigval:
                    min_eigval = eigval
                if eigval > max_eigval:
                    max_eigval = eigval
            
        l_normalizedeigenvalues = []
        for i in range(len(datadic['eigenvalues_perturbed']['data'][mode])):
            eigenvalue = sum(datadic['eigenvalues_perturbed']['data'][mode][i])/(len(datadic['eigenvalues_perturbed']['data'][mode][i])-1)
            eigenvalue_normalized = (eigenvalue-min_eigval)/(max_eigval-min_eigval)
            l_normalizedeigenvalues += [eigenvalue_normalized]

        plotname = 'average change in overlap'
        title1 = 'mode %s' %(mode+1)
        title = '"%s, %s \\n job: %s, chains: %s, cutoff: %s{\\305}"' %(plotname, title1, jobid, chains, cutoff_distance)
        gnuplot_plot(
            jobid, cutoff_distance, chains, title1,
            xtitle = 'residue', ytitle = 'average change in overlap', plotname = plotname,
            filename = '%s_deltaoverlap_mode%s' %(jobid, str(mode+1).zfill(2)),
            data = l_overlaps,
            title = title,
            )
        ## write overlaps to txt file
        lines = []
        for i in range(len(l_overlaps)):
            overlap = l_overlaps[i]
            lines += ['%3i %f\n' %(i, overlap)]
        fd = open('delta_overlaps_mode%s.txt' %(str(mode+1).zfill(2)),'w')
        fd.writelines(lines)
        fd.close()
        ## write overlaps to pdb file
        write_overlaps_to_pdb_file(
            d_coordinates, d_hessian, l_overlaps, mode, jobid, winsize = winsize,
            )
        
        ## write eigenvalue to pdb file
        write_overlaps_to_pdb_file(
            d_coordinates, d_hessian, l_normalizedeigenvalues, mode, jobid, winsize = winsize,
            s_data = 'eigenvalues',
            )

        ## group data
        d_overlaps = {
            'ALA':[],'CYS':[],'ASP':[],'GLU':[],'PHE':[],'GLY':[],'HIS':[],
            'ILE':[],'LYS':[],'LEU':[],'MET':[],'ASN':[],'PRO':[],'GLN':[],
            'ARG':[],'SER':[],'THR':[],'VAL':[],'TRP':[],'TYR':[],
            }
        if len(d_hessian.keys()) != len(l_overlaps):
            stop
        for i in range(len(l_overlaps)):
            res_name = d_hessian[i]['res_name']
            overlap = l_overlaps[i]
            d_overlaps[res_name] += [overlap]

        fileprefix = '%s_deltaoverlap_residue_mode%s' %(jobid, str(mode+1).zfill(2)),
        gnuplot_histogram(d_overlaps, fileprefix)

    axistitle = 'residue'
    ## set the diagonal elements to a standard value to get full contour
    for key in datadic.keys():
        for mode in range(6,12):
            if datadic[key]['diagonal'] == 'max':
                diagonal = max(datadic[key]['data'][mode][(winsize-1)/2])
            elif datadic[key]['diagonal'] == 'mode':
                diagonal = mode+1
            else:
                diagonal = datadic[key]['diagonal']
            for i in range(0+(winsize-1)/2,N-(winsize-1)/2):
                datadic[key]['data'][mode][i][i] = diagonal
                for j in range(
                    max(0+(winsize-1)/2,i-winsize+1),min(N-(winsize)/2,i+winsize-1)
                    ):
                    datadic[key]['data'][mode][i][j] = diagonal
##                for j in range(max(0,i-(winsize-1)/2),min(N-1,i+(winsize-1)/2)+1):
##                    datadic[key]['data'][mode][i][j] = diagonal
##            stop

##    ##
##    ## set titles
##    ##
##    xtitle = axistitle
##    ytitle = axistitle
##    ztitle = '{/Symbol D}overlap'
##    fileprefix = jobid+'overlapscombined_mode3N'
##    title = '%s, %s\\njob: %s\\nchains: %s\\ndistance cutoff: %s{\305}' %('{/Symbol D}overlap', 'modes 7-3N', jobid, chains, cutoff_distance),
##
##    ## overlaps_combined between nonperturbed and perturbed modes 7-3N
##    gnuplot_splot(
##        xtitle, ytitle, ztitle,
##        fileprefix,
##        datadic['overlaps_combined']['data'][-1],
##        title,
##        z1 = '', z2 = '',
##        )
##
##    ## add topology to margin of cross-correlations maps
##    fileprefix = jobid+'overlapscombined_mode3N'
##    topology(
##        fileprefix,chains,d_secondary,N,
##        )
    
    for key in datadic:
        print key
        for mode in range(6,12):
            xtitle = axistitle,
            ytitle = axistitle,
            ztitle = datadic[key]['name']
            fileprefix = '%s_%s_mode%s' %(jobid, key, str(mode+1).zfill(2))
            title = '%s, %s \\n job: %s, chains: %s, cutoff: %s{\\305}' %(datadic[key]['name'], mode+1, jobid, chains, cutoff_distance),
            gnuplot_splot(
                xtitle, ytitle, ztitle,
                fileprefix,datadic[key]['data'][mode],title,
                z1 = datadic[key]['zrange'][0], z2 = datadic[key]['zrange'][1],
                )

            ## add topology to margin of cross-correlations maps
            topology(
                fileprefix,chains,d_secondary, N,
                )

    return


def write_overlaps_to_pdb_file(
    d_coordinates, d_hessian, l_overlaps, mode, job, winsize = 1,
    s_data = 'overlaps',
    ):

    l_is = d_hessian.keys()
    l_is.sort()
    lines = []
    for i in l_is:
        chain = d_hessian[i]['chain']
        res_no = d_hessian[i]['res_no']
        iCode = d_hessian[i]['iCode']
        altloc = d_hessian[i]['altloc']
        res_name = d_hessian[i]['res_name']

        atom_names = d_coordinates['chains'][chain]['residues'][res_no]['iCodes'][iCode]['res_names'][res_name]['atoms'].keys()
        atom_names.sort()
        for atom_name in atom_names:
            altloc = min(d_coordinates['chains'][chain]['residues'][res_no]['iCodes'][iCode]['res_names'][res_name]['atoms'][atom_name]['altlocs'].keys())
            coordinate = d_coordinates['chains'][chain]['residues'][res_no]['iCodes'][iCode]['res_names'][res_name]['atoms'][atom_name]['altlocs'][altloc]['coordinate']
            atom_no = d_coordinates['chains'][chain]['residues'][res_no]['iCodes'][iCode]['res_names'][res_name]['atoms'][atom_name]['altlocs'][altloc]['atom_no']

            x = coordinate[0]
            y = coordinate[1]
            z = coordinate[2]


            ## append atom line
            record = 'ATOM'
            occupancy = 1.0
            if i < 0+(winsize-1)/2:
                bfactor = 100-100*l_overlaps[(winsize-1)/2]
            elif i > max(l_is)-(winsize-1)/2:
                bfactor = 100-100*l_overlaps[max(l_is)-(winsize-1)/2]
##            if i < 0+winsize+2:
##                bfactor = 100*l_overlaps[winsize+2]
##            elif i > max(l_is)-winsize-2:
##                bfactor = 100*l_overlaps[max(l_is)-winsize-2]
            else:
                bfactor = 100*l_overlaps[i]
            line = '%6s%5i  %3s%s%3s %1s%4i%1s   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s  \n' %(
                record.ljust(6),atom_no,
                atom_name.ljust(3), altloc, res_name, chain, res_no, iCode,
                x,y,z, occupancy, bfactor,atom_name[0].rjust(2)
                )
            lines += [line]

    fd = open('%s_%s_%s.pdb' %(job, s_data, str(mode+1).zfill(2)), 'w')
    fd.writelines(lines)
    fd.close()

    return


def pre_perturbation_plot(
    jobid,
    cutoff_distance,
    chains,d_secondary,
    eigenvectors_nonperturbed,
    eigenvectors_combined_nonperturbed,
    ):

    import math, Numeric

    for plottype in ['individual','combined',]:

        if plottype == 'individual':
            eigenvectors = eigenvectors_nonperturbed
        elif plottype == 'combined':
            eigenvectors = eigenvectors_combined_nonperturbed

        for mode in range(6,12):

            print 'plotting mode %s, %s' %(mode+1,plottype)

            ## calculate residue fluctuation (atom,RMSF) for each mode

            lengths = []
            for i in range(0,len(eigenvectors[mode]),3):
                length = math.sqrt(
                    eigenvectors[mode][i+0]**2+
                    eigenvectors[mode][i+1]**2+
                    eigenvectors[mode][i+2]**2
                    )
                lengths.append(length)

            ## write lengths to file

            lines = []
            for i in range(len(lengths)):
                lines.append(str(lengths[i])+'\n')
            fd = open('%s_%s_%s_lengths.txt' %(jobid,plottype,mode+1), 'w')
            fd.writelines(lines)
            fd.close()

            ## plot residue fluctuation (atom,RMSF) for each mode

            if plottype == 'individual':
                title1 = 'mode %s' %(mode+1)
            elif plottype == 'combined':
                title1 = 'modes 7-%s' %(mode+1)

            plotname = 'residue displacements'
            title = '"%s, %s \\n job: %s, chains: %s, cutoff: %s{\\305}"' %(plotname, title1, jobid, chains, cutoff_distance)
            gnuplot_plot(
                jobid, cutoff_distance, chains, title1,
                xtitle = 'residue', ytitle = 'RMSF', plotname = plotname,
                filename = '%s_%sdisplacement_mode%s' %(jobid, plottype, str(mode+1).zfill(2)),
                data = lengths, title = title,
                )

            ## prepare data for cross correlation maps
            ccdata = Numeric.zeros((len(eigenvectors[mode])/3,len(eigenvectors[mode])/3), typecode='d')
            for ccrow in range(len(ccdata)):
                for cccol in range(len(ccdata)):
                    ccdata[ccrow][cccol] = cosangle(eigenvectors[mode][3*ccrow:3*ccrow+3], eigenvectors[mode][3*cccol:3*cccol+3])

            ## plot cross-correlation map (residue,residue,cross-correlation) for each mode
            fileprefix = '%s_%s_crosscorrelation_mode%s' %(jobid, plottype, str(mode+1).zfill(2))

            for i in range(len(ccdata)):
                ccdata[i][i] = 1
            xtitle = 'residue'
            ytitle = 'residue'
            ztitle = 'cross correlation'
            title = '%s, %s \\n job: %s, chains: %s, cutoff: %s{\\305}' %('cross correlation map', mode+1, jobid, chains, cutoff_distance),

            gnuplot_splot(
                xtitle, ytitle, ztitle,
                fileprefix,ccdata,title,
                z1 = -1, z2 = 1,
                )

            ## add topology to margin of cross-correlations maps
            topology(
                fileprefix, chains, d_secondary, len(ccdata),
                )

    return


def gnuplot_splot(
    xtitle, ytitle, ztitle, fileprefix, data,
    title,
    z1 = None, z2 = None,
    ):

    import os

    ## write gnuplot data to txt file
    gnuplot_splot_data = []
    for x in range(len(data)):
        for y in range(len(data[x])):
            if x >= y:
                z = data[x][y]
            if y > x:
                z = data[y][x]
            gnuplot_splot_data.append('%4i %4i %16.13f\n' %(x+1, y+1, z))
        gnuplot_splot_data.append('%4i %4i %16.13f\n\n' %(x+1, len(data[x])+1, data[0][0]))
    for i in range(len(data)+1):
        gnuplot_splot_data.append('%4i %4i %16.13f\n' %(len(data)+1, i+1, data[0][0]))
    fd = open('%s.gnuplotdata' %(fileprefix), 'w')
    fd.writelines(gnuplot_splot_data)
    fd.close()
    ## write gnuplot settings to txt file
    lines = ['set size square\n'] ## scale square
    lines += [
        'set terminal postscript eps enhanced color "Helvetica" 48\n',
        'set output "%s.ps"\n' %(fileprefix),
        'set size 4,4\n', ## scale 400%
        'set view map\n', ## change orientation of plot
        'set autoscale fix\n', ## scale axes
        'set style data pm3d\n', ## set by default?
        'set style function pm3d\n', ## set by default?
        'set encoding iso_8859_1\n', ## postscript encoding for special characters
##        'set title %s, 4\n' %(title),
        'set title "%s" offset 0,1\n' %(title),
        'set xlabel "%s"\n' %(xtitle),
        'set ylabel "%s"\n' %(ytitle),
        'set cblabel "%s"\n' %(ztitle),
        ]
    if z1 != None and z2 != None:
        lines += [
            'set cbrange [%s:%s]\n' %(z1,z2), ## set colorbox range
            ]
    lines += [
        'set palette model CMY rgbformulae 7,5,15\n',
        'set pm3d map corners2color c1\n', ## generate a 2D surface rather than 3D points
        'splot "%s.gnuplotdata" title ""\n' %(fileprefix), ## splot gnuplot data file
        ## set ticslevel
        ]

    fd = open('%s.gnuplotsettings' %(fileprefix), 'w')
    fd.writelines(lines)
    fd.close()
    ## plot data with gnuplot splot
    os.system('/software/bin/gnuplot %s.gnuplotsettings' %(fileprefix))
    ## convert postscript to portable network graphics
    os.system('convert %s.ps %s.png' %(fileprefix, fileprefix))
    os.remove('%s.ps' %(fileprefix))
    os.remove('%s.gnuplotdata' %(fileprefix))
    os.remove('%s.gnuplotsettings' %(fileprefix))

    return


def gnuplot_plot(jobid, cutoff_distance, chains, title1, xtitle, ytitle, plotname, filename, data, title):

    import os

    ## write gnuplot data to txt file
    lines = []
    lines.append('%4i %16.13f\n' %(0, 0))
    for i in range(len(data)):
        lines.append('%4i %16.13f\n' %(i+1, data[i]))
    lines.append('%4i %16.13f\n' %(len(data)+1, 0))
    fd = open('%s.gnuplotdata' %(filename), 'w')
    fd.writelines(lines)
    fd.close()
    ## write gnuplot settings to txt file
    lines = [
        'set terminal postscript eps enhanced color "Helvetica" 48\n',
        'set output "%s.ps"\n' %(filename),
        'set size 4,4\n', ## scale 400%
        'set autoscale fix\n', ## scale axes
        'set encoding iso_8859_1\n', ## postscript encoding for special characters
        'set style data histeps\n', ## change plot style to histogram
        'set style function histeps\n', ## change plot style to histogram
##        'set title %s ,4\n' %(title),
        'set title %s offset 0,1\n' %(title),
        'set xlabel "%s"\n' %(xtitle),
        'set ylabel "%s"\n' %(ytitle),
        'plot "%s.gnuplotdata" title "" lt 3 lw 16\n' %(filename), ## plot gnuplot data file
        ## set ticslevel
        ]

    fd = open('%s.gnuplotsettings' %(filename), 'w')
    fd.writelines(lines)
    fd.close()
    ## plot data with gnuplot plot
    os.system('/software/bin/gnuplot %s.gnuplotsettings' %(filename))
    ## convert postscript to portable network graphics
    os.system('convert %s.ps %s.png' %(filename, filename))
    os.remove('%s.ps' %(filename))
    os.remove('%s.gnuplotdata' %(filename))
    os.remove('%s.gnuplotsettings' %(filename))

    return


def calculate_distance_matrix(l_coordinates):

    import Numeric

    N = len(l_coordinates)

    matrix_distances = Numeric.zeros((N,N),typecode='d')
    for row in range(N):
        for col in range(N):
            if row > col:
                v = l_coordinates[row]-l_coordinates[col]
                matrix_distances[row][col] = sum(v**2)
                matrix_distances[col][row] = sum(v**2)

    return matrix_distances


def parse_dictionary_of_coordinates(d_coordinates, chains, atoms_hessian):

    '''this function was previously named calculate_N()'''

    i = 0
    l_coordinates = []
    d_hessian = {}

    for chain in chains:
        ## assume sequential numbering of residues
        res_nos = d_coordinates['chains'][chain]['residues'].keys()
        res_nos.sort()
        for res_no in res_nos:
            ## assume sequential iCodes
            iCodes = d_coordinates['chains'][chain]['residues'][res_no]['iCodes'].keys()
            iCodes.sort()
            for iCode in iCodes:
                altloc_residue = min(d_coordinates['chains'][chain]['residues'][res_no]['iCodes'][iCode]['altlocs'].keys())
                res_name = d_coordinates['chains'][chain]['residues'][res_no]['iCodes'][iCode]['altlocs'][altloc_residue]['res_name']
                atom_names = d_coordinates['chains'][chain]['residues'][res_no]['iCodes'][iCode]['res_names'][res_name]['atoms'].keys()
                atom_names.sort()
                for atom_name in atom_names:
                    ## assume minimum/first altloc to be best alternative (highest occupancy)
                    altloc = min(d_coordinates['chains'][chain]['residues'][res_no]['iCodes'][iCode]['res_names'][res_name]['atoms'][atom_name]['altlocs'].keys())
                    coordinate = d_coordinates['chains'][chain]['residues'][res_no]['iCodes'][iCode]['res_names'][res_name]['atoms'][atom_name]['altlocs'][altloc]['coordinate']
                    record = d_coordinates['chains'][chain]['residues'][res_no]['iCodes'][iCode]['res_names'][res_name]['atoms'][atom_name]['altlocs'][altloc]['record']
                    atom_no = d_coordinates['chains'][chain]['residues'][res_no]['iCodes'][iCode]['res_names'][res_name]['atoms'][atom_name]['altlocs'][altloc]['atom_no']
                    element = d_coordinates['chains'][chain]['residues'][res_no]['iCodes'][iCode]['res_names'][res_name]['atoms'][atom_name]['altlocs'][altloc]['element']
                    if atom_name in atoms_hessian:

                        d_coordinates['chains'][chain]['residues'][res_no]['iCodes'][iCode]['res_names'][res_name]['atoms'][atom_name]['altlocs'][altloc]['i'] = i
                        ## append to d_hessian
                        d_hessian[i] = {
                            'chain':chain,'res_no':res_no,'iCode':iCode,'atom_name':atom_name,
                            'altloc':altloc,'res_name': res_name,'coordinate':coordinate,
                            'record':record,'atom_no':atom_no,'element':element,
                            }
                        ## append to l_coordinates
                        l_coordinates += [coordinate]

                        i += 1

    N = len(l_coordinates)

    return N, d_hessian, l_coordinates


def fill_matrix_hessian(matrix_hessian, factor, row_sup, col_sup, dist_sq, vector):

    for row_sub in range(3):
        for col_sub in range(3):

            if col_sub >= row_sub: #fill supersymmetrical elements when j>=i
                if dist_sq == 0:
                    print row_sub, col_sub, row_sup, col_sup
                value = factor*-vector[row_sub]*vector[col_sub]/dist_sq
                matrix_hessian[3*row_sup+row_sub,3*col_sup+col_sub] = value ##upper super off-diagonal; xixj, xiyj, xizj, yiyj, yizj, zizj
                matrix_hessian[3*col_sup+col_sub,3*row_sup+row_sub] = value ##lower super off-diagonal; xjxi, yjxi, zjxi, yjyi, zjyi, zjzi
                matrix_hessian[3*row_sup+row_sub,3*row_sup+col_sub] -= value ##super diagonal (row); xixi, xiyi, xizi, yiyi, yizi, zizi
                matrix_hessian[3*col_sup+col_sub,3*col_sup+row_sub] -= value ##super diagonal (col); xjxj, yjxj, zjxj, yjyj, zjyj, zjzj
                if col_sub > row_sub: #fill lower subsymmetrical elements
                    matrix_hessian[3*row_sup+col_sub,3*col_sup+row_sub] = value #upper super off-diagonal; yixj, zixj, ziyj
                    matrix_hessian[3*col_sup+row_sub,3*row_sup+col_sub] = value #lower super off-diagonal; xjyi, xjzi, yjzi
                    matrix_hessian[3*row_sup+col_sub,3*row_sup+row_sub] -= value ##super diagonal; yixi, zixi, ziyi
                    matrix_hessian[3*col_sup+row_sub,3*col_sup+col_sub] -= value ##super diagonal; yjxj, zjxj, zjyj

    return matrix_hessian

if __name__=='__main__':
    None
