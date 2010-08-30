#!/usr/bin/env python
##
##$Id$
##
##Tommy Carstensen, University College Dublin, 2005

import sys, Utilities
sys.path.append("../")
import Environment

class vibration:

    def main(
        self, jobid, filename,
        atoms_hessian = ['CA'], frames = 50,
        cutoff_distance = 10.,
        path_python = None, verbose = False, paralleldir = '',
        biomolecule = None, chains = [], model = None,
        residueIndexes = [],
        pre_perturbation_plot = True,
        ):

        '''
        Use first model if no model specified by user.
        chain(s): Y, biomolecule: Y; parse chains specified by user and apply transformation
        chain(s): Y, biomolecule: N; parse chains specified by user but don't apply transformation
        chain(s): N, biomolecule: Y; parse chains of biomolecule and apply transformation
        chain(s): N, biomolecule: N; parse chains of first biomolecule and apply transformation
        '''

        import os, Numeric, Utilities

        results = []

        fd = open(filename,'r')
        lines = fd.readlines()
        fd.close()

        ## parse pdb
        (
            d_REMARK350,
            d_primary, ## i.e. SEQRES, MODRES
            d_secondary, ## i.e. HELIX, SHEET
            d_coordinates, ## i.e. ATOM, HETATM, TER, MODEL, ENDMDL
            d_ligands,
            ) = Utilities.parse_pdb(lines, chains)

        ## assume multimeric biological unit if chains not specified by user
        if chains == []:
            chains = d_coordinates['chains'].keys()
            chains.sort()

        ##
        ## calculate N and convert coordinates from dic to list
        ##
        N, d_hessian, l_coordinates = Utilities.parse_dictionary_of_coordinates(d_coordinates, chains, atoms_hessian)

        ##
        ## calculate distance matrix
        ##
        matrix_distances = Utilities.calculate_distance_matrix(l_coordinates)

        ##
        ## calculate hessian matrix
        ##
        matrix_hessian = self.hessian_calculation(N, d_coordinates, chains, atoms_hessian, float(cutoff_distance), d_secondary, matrix_distances, l_coordinates, verbose = verbose)

        ##
        ## diagonalize hessian matrix
        ##
        eigenvectors_nonperturbed, eigenvalues_nonperturbed, eigenvectors_comb_nonperturbed, = Utilities.eigenv_calccomb(
                matrix_hessian, jobid, verbose,
                )

        ##
        ## set data lists and append matrices to be plotted for each combination of modes 6-12 before initiating loops over the two axes of the plot
        ##
        datadic = Utilities.datadic_return(N)

        ##
        ## loop over remres1
        ##

	environment = Environment.Environment()
	array = environment.splitArray(residueIndexes)
	winsize = 1
        for remres1 in range(N):

            ##
            ## loop over remres2
            ##
            for remres2 in range(N):

		if remres2 <= remres1:
			continue

		if not (remres1 in array or remres2 in array):
                    continue

                environment.output("[MODES] Residues: %s %s" % (remres1, remres2), rootOnly=False)

                l_rem = []
                for i in range(-(winsize-1)/2,(winsize+1)/2):
                    l_rem += [remres1+i,remres2+i]
                l_rem.sort()

                ##
                ## remove selected alpha carbon atoms!!!
                ##
                l_coordinates_perturbed = l_coordinates[:remres1-(winsize-1)/2]+l_coordinates[remres1+1+(winsize-1)/2:remres2-(winsize-1)/2]+l_coordinates[remres2+1+(winsize-1)/2:]
            
                matrix_hessian_perturbed = self.hessian_calculation(
                    N-len(l_rem), d_coordinates, chains, atoms_hessian, float(cutoff_distance), d_secondary, matrix_distances, l_coordinates_perturbed, verbose = verbose,
                    )

                (
                    eigenvectors_perturbed, eigenvalues_perturbed, eigenvectors_perturbed_combined,
                    ) = Utilities.eigenv_calccomb(
                        matrix_hessian_perturbed, jobid, verbose
                        )

                (
                    overlaps_single, max_overlaps_single, perturbed_modes_of_max_overlap_single, delta_perturbed_eigenvalues_of_max_overlap,
                    ) = Utilities.overlap_calculation(
                        eigenvectors_perturbed, eigenvectors_nonperturbed, eigenvalues_perturbed, eigenvalues_nonperturbed, l_rem
                        )

                (
                    overlaps_combined, max_overlaps_combined, perturbed_modes_of_max_overlap_combined,
                    ) = Utilities.overlap_calculation(
                        eigenvectors_perturbed_combined, eigenvectors_nonperturbed, eigenvalues_perturbed, eigenvalues_nonperturbed, l_rem,
                        )[:-1]

                datadic_loop = {
                    'eigenvalues_perturbed': eigenvalues_perturbed,
                    'emo': delta_perturbed_eigenvalues_of_max_overlap,
                    'overlaps_single': overlaps_single,
                    'overlaps_max': max_overlaps_single,
                    'mmo': perturbed_modes_of_max_overlap_single,
                    'overlaps_combined': overlaps_combined
                    }
                for key in datadic:
                    for mode in range(6,12):
                        datadic[key]['data'][mode][remres1][remres2] = datadic_loop[key][mode]
                        datadic[key]['data'][mode][remres2][remres1] = datadic_loop[key][mode]
                datadic['overlaps_combined']['data'][-1][remres1][remres2] = overlaps_combined[-1]
                datadic['overlaps_combined']['data'][-1][remres2][remres1] = overlaps_combined[-1]

                environment.output("[MODES] Overlap %lf" % (overlaps_single[6]), rootOnly=False)

        l_averages = []
        for remres1 in range(len(data =  datadic['overlaps_single']['data'][6])):
            l_overlaps = []
            for remres2 in range(len(datadic['overlaps_single']['data'][6][remres1])):
                overlap = datadic['overlaps_single']['data'][6][remres1][remres2]
                l_overlaps += [overlap]
            l_averages += [sum(l_overlaps)/len(l_overlaps)]

	print 'I am processors %d - results %s' % (environment.rank(), l_averages)

	results = environment.combineArray(l_averages)

	if environment.isRoot():
		print 'Combined results %s' % results

        return results


    def hessian_calculation(self, N, d_coordinates, chains, atoms_hessian, cutoff, d_secondary, matrix_distances, l_coordinates, l_rem = [], verbose = False):

        '''This function calculates the Hessian matrix. At least its supposed to...'''

        if verbose == True:
            print 'calculating the hessian/second-order derivate Kirchoff/Laplacian matrix'
        
        import math, Numeric

        cutoff_sq = cutoff**2

        matrix_hessian = Numeric.zeros((3*N,3*N), typecode='d')
        
        for row_sup in range(N):
            for col_sup in range(N):

                if col_sup > row_sup:
                    xi = l_coordinates[row_sup][0]
                    xj = l_coordinates[col_sup][0]
                    yi = l_coordinates[row_sup][1]
                    yj = l_coordinates[col_sup][1]
                    zi = l_coordinates[row_sup][2]
                    zj = l_coordinates[col_sup][2]
                    x = xj-xi
                    y = yj-yi
                    z = zj-zi
                    dist_sq = x**2+y**2+z**2
##                    if dist_sq <= cutoff_sq:
                    sigmoidfactor = Utilities.sigmoid(math.sqrt(dist_sq), cutoff)
                    vector = [x,y,z]
                    for row_sub in range(3):
                        for col_sub in range(3):

                            if col_sub >= row_sub: #fill supersymmetrical elements when j>=i
                                value = sigmoidfactor*-vector[row_sub]*vector[col_sub]/dist_sq
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
        
if __name__ == '__main__':
    filename = sys.argv[-1]
    jobid = 'PEATSA'
    instance_vibration = vibration()
    instance_vibration.main(jobid,filename,residueIndexes=range(0,4))
