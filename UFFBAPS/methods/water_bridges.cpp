/*
 #
 # UFFBAPS - Unified Force Field for Binding And Protein Stability
 # Copyright (C) 2010 Jens Erik Nielsen & Chresten Soendergaard
 #
 # This program is free software: you can redistribute it and/or modify
 # it under the terms of the GNU General Public License as published by
 # the Free Software Foundation, either version 3 of the License, or
 # (at your option) any later version.
 #
 # This program is distributed in the hope that it will be useful,
 # but WITHOUT ANY WARRANTY; without even the implied warranty of
 # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 # GNU General Public License for more details.
 #
 # You should have received a copy of the GNU General Public License
 # along with this program.  If not, see <http://www.gnu.org/licenses/>.
 #
 # Contact information: 
 # Email: Jens.Nielsen_at_gmail.com
 # Normal mail:
 # Jens Nielsen
 # SBBS, Conway Institute
 # University College Dublin
 # Dublin 4, Ireland
 */


#include "water_bridges.h"


Water_Bridges::Water_Bridges():Method()
{
  set_parameters();  
}

Water_Bridges::Water_Bridges(FILE * reporter):Method(reporter)
{
  /*** write out to out file **/
  writeDescription(reporter);
  set_parameters();  
}

Water_Bridges::Water_Bridges(FILE * reporter, FILE * resultsFile, vector<Soup *> A, vector<Soup *> B):Method(reporter)
{
  /*** write out to out file **/
  writeDescription(reporter);
  set_parameters();  

  double result = 0;

  if(A.size() != B.size())
    {
      cout << "WARNING: An equal amount of protein and ligand structures must be given to the water bridge function!\n";
      return;
    }
  
  for(unsigned int i=0; i<A.size(); i++)
    {
      result = calculate(A[i],B[i]);
      printf("Water bridge for %s-%s gave %8.2f\n",A[i]->name.c_str(), B[i]->name.c_str(), result);
      
      string decoy_number;
      if(B[i]->name.find_last_of("_")== string::npos)
      	decoy_number = "-1";
      else
	decoy_number = B[i]->name.substr(B[i]->name.find_last_of("_")+1);
      
       fprintf(resultsFile,"%-4s %-4s %8f\n",
	      A[i]->name.substr(0,4).c_str(),
	      decoy_number.c_str(),
	      result);
      
    }  

}

Water_Bridges::Water_Bridges(FILE * resultsFile, FILE * reporter, vector<Soup *> A):Method(reporter)
{
  /*** write out to out file **/
  writeDescription(reporter);
  set_parameters();  
    

  double result = 0;

  for(unsigned int i=0; i<A.size(); i++)
    {
      result = calculate(A[i]);
      printf("Water bridge for %s gave %8.2f\n",A[i]->name.c_str(), result);
    }  

}


Water_Bridges::~Water_Bridges(){}




float Water_Bridges::calculate(Soup * A)
{
  potential_waters.clear();

  if(!A->bonds_assigned)
    BM.assign_bonds(A);

  // place potential water oxygen atoms
  dry_atoms = A->get_dry_atoms();
  place_water_contacts();

  // count number of water bridges
  float res = 0.0;
  vector<potential_water_struct>::iterator water;
  vector<Atom*>::iterator contact1, contact2;

  //Atom * c1,*c2;

  int max_res_dist,res_dist;
  float this_contribution, sum_desolv_factors;
  // for each placed water
  for(water  = potential_waters.begin(); water != potential_waters.end(); water++)
    if(water->active)
      if (water->contacts.size() > 1)
	{
	  //cout<<water->contacts.size()<<" contacts for ";///////////////////////
	  //water->water_atom.display();////////////////////////////

	  sum_desolv_factors = 0;
	  // for each contact
	  for(contact1  = water->contacts.begin(); contact1 != water->contacts.end(); contact1++)
	    {
	      this_contribution = water_bridge_energy;
	      max_res_dist = 0;
	      
	      // find biggest residue dist to other contacts
	      for(contact2  = water->contacts.begin(); contact2 != water->contacts.end(); contact2++)
		{
		  res_dist = abs((*contact1)->residue_no - (*contact2)->residue_no);
		  if (max_res_dist<res_dist)
		    max_res_dist = res_dist;
		}
	      
	      // find the contribution of this contact to the water molecule
	      if (max_res_dist == 0)
		this_contribution = this_contribution * same_residue_factor;
	      else if (max_res_dist == 1)
		this_contribution = this_contribution * neighbour_residue_factor;
	      
	      if(include_surf_factor == 1)
		this_contribution = this_contribution * (*contact1)->desolvation_factor;
	      
	      
	      //cout<<"   "<<this_contribution<<":  "<<(*contact1)->desolvation_factor<<endl;//////////////////////////////
	      //(*contact1)->display();/////////////////////////////
	      
	      res += this_contribution;
	      
	      sum_desolv_factors+= (*contact1)->desolvation_factor;
	    }
	  
	  
	  // add desolvation penalty for the water
	  res +=  water_desolvation_penalty * sum_desolv_factors/water->contacts.size();
      }

  //write_water_atoms(A->name);
  return res;
}



float Water_Bridges::calculate(Soup * A, Soup * B)
{

  potential_waters.clear();

  if(!A->bonds_assigned)
    BM.assign_bonds(A);

  if(!B->bonds_assigned)
    BM.assign_bonds(B);

  // remember to get waters from both soups
  dry_atoms_A = A->get_dry_atoms();
  dry_atoms_B = B->get_dry_atoms();

  dry_atoms = dry_atoms_A;
  dry_atoms.insert(dry_atoms.end(), dry_atoms_B.begin(), dry_atoms_B.end());


  // place potential water oxygen atoms
  place_water_contacts();

  float res = 0.0;
  vector<potential_water_struct>::iterator water;
  vector<Atom*>::iterator contact;
  float this_contribution, surf_factor;
  int no_protein_contacts, no_ligand_contacts;
  // for each placed water

  for(water  = potential_waters.begin(); water != potential_waters.end(); water++)
    if(water->active)
      {
	// try to find waters with contacts in both protein and ligand
	this_contribution = 0.0;
	no_protein_contacts = 0;
	no_ligand_contacts = 0;
	surf_factor=0;
	for(contact = water->contacts.begin(); contact != water->contacts.end(); contact++)
	  {
	    if( find(dry_atoms_A.begin(), dry_atoms_A.end(), *contact) != dry_atoms_A.end())
	      no_protein_contacts++;
	    if( find(dry_atoms_B.begin(), dry_atoms_B.end(), *contact) != dry_atoms_B.end())
	      no_ligand_contacts++;
	    
	    surf_factor += (*contact)->desolvation_factor;
	    
	  }
	
	if (no_protein_contacts > 0 && no_ligand_contacts > 0)//this water has contact in both protein and water
	  {
	    this_contribution = water_bridge_energy*(no_protein_contacts+no_ligand_contacts); // contribution is number of contacts
	    if(include_surf_factor == 1)
	      this_contribution = water_bridge_energy*surf_factor; //contribution is sum of surf_factors

	    res += this_contribution;

	    
	    // add desolvation penalty for the water
	    res +=  water_desolvation_penalty * surf_factor/((float) water->contacts.size());
	    

	    //water->water_atom.display();
	    //cout<<"Contributes "<<this_contribution<<" and the water has desolvation penalty "<<water_desolvation_penalty * surf_factor/water->contacts.size()
	    //<<" and "<<water->contacts.size()<<" contacts"<< endl;

	  }
	//	else
	//  water->active = false;
      }


  //write_water_atoms(A->name);
  return res;

}
  

void Water_Bridges::find_water_contacts(Soup * A)
// this will mark all acceptor/donor atoms able to contact a water molecule
{ 
  potential_waters.clear();

  if(!A->bonds_assigned)
    BM.assign_bonds(A);

  // place potential water oxygen atoms
  dry_atoms = A->get_dry_atoms();
  place_water_contacts();


  // mark atoms with water contact
  vector<potential_water_struct>::iterator water;
  vector<Atom*>::iterator contact;

  // for each placed water
  for(water  = potential_waters.begin(); water != potential_waters.end(); water++)
      for(contact = water->contacts.begin(); contact != water->contacts.end(); contact++)
	(*contact)->hydrogen_bond_could_be_replaced_by_water = true;

  return;
}



void Water_Bridges::place_water_contacts()
{

  potential_water_boxes.clear();
  potential_water_boxes_incl_neighbours.clear();

  // place all possible water atoms
  vector<Atom*>::iterator atom;
  for(atom = dry_atoms.begin(); atom != dry_atoms.end(); atom++)
    {

      if(HB.is_potential_acceptor(*atom) or HB.is_potential_donor(*atom))
	if ( !no_heavy_atoms_bonded((*atom)) ) //to avoid e.g. free oxygens (in e.g. 1tni)
	  place_potential_water_contacts(*atom);
    }

  // divide potential waters into boxes
  make_boxes();
  //write_water_atoms("newv_ALL");
 
  // remove those too close to other atoms
  for(unsigned int i=0; i<potential_water_boxes.size();i++)
    remove_overlaps_with_protein(potential_water_boxes[i], atom_boxes_incl_neighbours[i]);

  //write_water_atoms("newv_NO_OVERLAP");


  // merge close potential water atoms
  for(unsigned int i=0; i<potential_water_boxes.size();i++)
    merge_waters(potential_water_boxes[i], potential_water_boxes_incl_neighbours[i]);

  //write_water_atoms("newv_MERGED");
  
  // remove water-water bumps
  for(unsigned int i=0; i<potential_water_boxes.size();i++)
    remove_water_bumps(potential_water_boxes[i], potential_water_boxes_incl_neighbours[i]);
  

  
  // remove those too close to other atoms
  //  for(unsigned int i=0; i<potential_water_boxes.size();i++)
  // remove_overlaps_with_protein(potential_water_boxes[i], atom_boxes_incl_neighbours[i]);


  // write_water_atoms("newv_waters");
  
  return;
}


void Water_Bridges::make_boxes()
{
  float 
    box_size = 5.0,
    minx= 10000.0,miny= 10000.0,minz= 10000.0,
    maxx=-10000.0,maxy=-10000.0,maxz=-10000.0;

  vector<potential_water_struct>::iterator potential_water;
  vector<Atom*>::iterator atom;

  // find min and max values
  for (potential_water=potential_waters.begin(); potential_water!=potential_waters.end(); potential_water++)
    {
      if(minx > potential_water->water_atom.x)
	minx = potential_water->water_atom.x;

      if(miny > potential_water->water_atom.y)
	miny = potential_water->water_atom.y;

      if(minz > potential_water->water_atom.z)
	minz = potential_water->water_atom.z;

      if(maxx < potential_water->water_atom.x)
	maxx = potential_water->water_atom.x;

      if(maxy < potential_water->water_atom.y)
	maxy = potential_water->water_atom.y;

      if(maxz < potential_water->water_atom.z)
	maxz = potential_water->water_atom.z;

    }

  // also check min and maxs for protein/ligand atoms
  for(atom = dry_atoms.begin(); atom != dry_atoms.end(); atom++)
    {
      if(minx > (*atom)->x)
	minx = (*atom)->x;

      if(miny > (*atom)->y)
	miny = (*atom)->y;

      if(minz > (*atom)->z)
	minz = (*atom)->z;

      if(maxx < (*atom)->x)
	maxx = (*atom)->x;

      if(maxy < (*atom)->y)
	maxy = (*atom)->y;

      if(maxz < (*atom)->z)
	maxz = (*atom)->z;

    }
  


  int 
    imaxx = (int) ( ceil((maxx - minx)/box_size)+1 ),
    imaxy = (int) ( ceil((maxy - miny)/box_size)+1 ),
    imaxz = (int) ( ceil((maxz - minz)/box_size)+1 ),
    this_ix,this_iy,this_iz;


  vector<int> temp;
  potential_water_boxes.assign(imaxx*imaxy*imaxz, temp);
  potential_water_boxes_incl_neighbours.assign(imaxx*imaxy*imaxz, temp);
  atom_boxes_incl_neighbours.assign(imaxx*imaxy*imaxz, temp);
  atom_boxes.assign(imaxx*imaxy*imaxz, temp);

  // put water atoms into boxes
  for (unsigned int i=0; i<potential_waters.size(); i++)
    {    
      this_ix = ((int) ( (potential_waters[i].water_atom.x - minx)/box_size ));
      this_iy = ((int) ( (potential_waters[i].water_atom.y - miny)/box_size ));
      this_iz = ((int) ( (potential_waters[i].water_atom.z - minz)/box_size ));

      potential_water_boxes.at(this_ix +
			       this_iy*imaxx +
			       this_iz*imaxx*imaxy).push_back(i);

    }

  // put dry atoms into boxes
  for (unsigned int i=0; i<dry_atoms.size(); i++)
    {    
      this_ix = ((int) ( (dry_atoms[i]->x - minx)/box_size ));
      this_iy = ((int) ( (dry_atoms[i]->y - miny)/box_size ));
      this_iz = ((int) ( (dry_atoms[i]->z - minz)/box_size ));

      atom_boxes.at(this_ix +
		    this_iy*imaxx +
		    this_iz*imaxx*imaxy).push_back(i);
      
    }

  // join neighbours
  int 
    startx,starty,startz,
    endx,endy,endz,
    index,this_index;

  for(int this_ix=0; this_ix<imaxx; this_ix++) // for each box
    for(int this_iy=0; this_iy<imaxy; this_iy++)
      for(int this_iz=0; this_iz<imaxz; this_iz++)
	{
	  startx = max(this_ix - 1, 0);
	  starty = max(this_iy - 1, 0);
	  startz = max(this_iz - 1, 0);
	  
	  endx = min(this_ix + 2, imaxx);
	  endy = min(this_iy + 2, imaxy);
	  endz = min(this_iz + 2, imaxz);

	  this_index = this_ix + this_iy*imaxx + this_iz*imaxx*imaxy;


	  for(int ix=startx; ix<endx; ix++)
	    for(int iy=starty; iy<endy; iy++)
	      for(int iz=startz; iz<endz; iz++)
		{
		  index = ix + iy*imaxx + iz*imaxx*imaxy;
		  (potential_water_boxes_incl_neighbours.at(this_index)).insert((potential_water_boxes_incl_neighbours.at(this_index)).end(),
										(potential_water_boxes.at(index)).begin(),
										(potential_water_boxes.at(index)).end());
		  
		  
		  (atom_boxes_incl_neighbours.at(this_index)).insert((atom_boxes_incl_neighbours.at(this_index)).end(),
								     (atom_boxes.at(index)).begin(),
								     (atom_boxes.at(index)).end());

		}
	}
}


void Water_Bridges::write_water_atoms(string name)
{

  // write out water atoms
  vector<Atom*> water_atoms;

  vector<potential_water_struct>::iterator potential_water;

  for(potential_water =  potential_waters.begin();
      potential_water != potential_waters.end();
      potential_water++)
    if(potential_water->active)
      {
	water_atoms.push_back(&(potential_water->water_atom));
	/*
	cout<<"Water ";
	potential_water->water_atom.display();
	//cout<<"  has "<<potential_water_contacts[i].contacts.size()<<" contacts"<<endl;
	for (int k=0; k<potential_water->contacts.size();k++)
	  {
	    cout<<"      ";
	    (potential_water->contacts[k])->display();
	  }
	*/
      }

  Write_File fw;
  fw.write_atoms(water_atoms, "pdb", "ATOM  ", name+"_water.pdb");



}



void Water_Bridges::place_potential_water_contacts(Atom * atom)
{
  //cout<<"Placing waters for: ";
  //atom->display();

  Vector root,pos,perp,norm;
  vector<Vector> water_vectors;

  root = find_root(atom);
  pos  = atom->get_pos();
  norm = pos-root;
 
  // get first vector to potential water
  perp = norm.get_perpendicular();

  norm.normalise();
  perp.normalise();

  // generate possible water vectors
  for(unsigned int i=0; i<d_along.size(); i++)
    water_vectors.push_back((norm*d_along[i]) + (perp*d_perp[i])+pos);

  // rotate
  float angle;

  for(unsigned int i=0; i< water_vectors.size(); i++)
    {
      
      // store original water vector
      angle = 0.0;
      potential_water_struct pws((Atom(water_vectors[i], atom->residue_no)));
      pws.add_contact(atom);

      potential_waters.push_back(pws);
      //all_waters.push_back(&(potential_waters[potential_waters.size()-1].water_atom));


				 
      while(angle<360)
	{

	  //cout<<"step:"<<steps[i]<<" at "<<i<<" in "<< water_vectors.size()<<endl;
	  //printf("%4.0f ",angle);
	  //water_vector->display("");

	  water_vectors[i].rotate(steps[i], root, pos);
	  if (fabs(water_vectors[i].x)>1000.0)
	    cout<<"FIS"<<endl;

	  // store rotated vectors
	  potential_water_struct pws((Atom(water_vectors[i], atom->residue_no)));
	  pws.add_contact(atom);
	        
	  potential_waters.push_back(pws);
	  //all_waters.push_back(&(potential_waters[potential_waters.size()-1].water_atom));

	  angle+=steps[i];
	}
    }


 
}




Vector Water_Bridges::find_root(Atom * atom)
{

  float count = 0;
  Vector res(0.0, 0.0, 0.0);

  vector<Atom*>::iterator bonded_atom;

  for(bonded_atom = atom->bonded_atoms.begin(); bonded_atom != atom->bonded_atoms.end(); bonded_atom++)
    if((*bonded_atom)->element != "H")
      {
	count = count + 1;
	res.x += (*bonded_atom)->x;
	res.y += (*bonded_atom)->y;
	res.z += (*bonded_atom)->z;
      }
  
  
  res.x = res.x/count;
  res.y = res.y/count;
  res.z = res.z/count;
  if (isnan(res.x))
    {
      printf("NAN in find root");
      atom->display();
	
    }


  return res;
}


void Water_Bridges::remove_overlaps_with_protein(vector<int> indices, vector<int> dry_neighbour_indices)
{
  vector<potential_water_struct>::iterator potential_water1, potential_water2;
  vector<Atom*>::iterator atom;

  // for each potential water 
  for(unsigned int i=0; i<indices.size();i++)
    if(potential_waters[indices[i]].active)
      {
	for(unsigned int d=0; d<dry_neighbour_indices.size(); d++)
	  {
	    if(!(dry_atoms[dry_neighbour_indices[d]]->is_hydrogen))
	    if(dist.calculate(dry_atoms[dry_neighbour_indices[d]], &(potential_waters[indices[i]].water_atom) , false) < protein_bump_dist_sq)
	      if(find(potential_waters[indices[i]].contacts.begin(),
		      potential_waters[indices[i]].contacts.end(),
		      dry_atoms[dry_neighbour_indices[d]] ) == potential_waters[indices[i]].contacts.end())
		{
		  //cout<<"removing "<<indices[i]<<endl;
		  potential_waters[indices[i]].active = false;
		  break; //no reason to look for more atoms
		  
		}
	  }
	
      }
	/*
  for(potential_water1  = potential_waters.begin(); 
      potential_water1 != potential_waters.end();)
    {
      overlap = false;
      for(atom = atoms.begin(); atom != atoms.end(); atom++)
	if(!(*atom)->is_hydrogen)
	  if (fabs((*atom)->x - potential_water1->water_atom.x) < protein_bump_dist)
	    if (fabs((*atom)->y - potential_water1->water_atom.y) < protein_bump_dist)
	      if (fabs((*atom)->z - potential_water1->water_atom.z) < protein_bump_dist)
		if(dist.calculate(*atom, &(potential_water1->water_atom) , false) < protein_bump_dist_sq)
		  if(find(potential_water1->contacts.begin(),
			  potential_water1->contacts.end(),
			  *atom) == potential_water1->contacts.end())
		    {
		      overlap = true;
		      break; //no reason to look for more atoms

		    }
      
      if (overlap)
	  potential_water1->active = false;

      potential_water1++;
      
    }
	*/
  return;
}




void Water_Bridges::merge_waters(vector<int> indices, vector<int> neighbour_indices)
{
   
  vector<Atom*>::iterator contact1,contact2;
  bool shared_contact;
  
  //  cout<<"Merging "<<indices.size()<<" waters"<<endl;

  for(unsigned int i=0; i<indices.size();i++)
    if(potential_waters[indices[i]].active)
      for(unsigned int j=0; j<neighbour_indices.size();j++)
	if(indices[i] != neighbour_indices[j])
	  if(potential_waters[neighbour_indices[j]].active)
	    {
	      shared_contact = false; //first check if waters have same primary contact
	      if(potential_waters[indices[i]].contacts[0] == potential_waters[neighbour_indices[j]].contacts[0])
		shared_contact = true;
	      
	      if(!shared_contact) //if not go through all contacts
		for(contact1  = potential_waters[indices[i]].contacts.begin();
		    contact1 != potential_waters[indices[i]].contacts.end();
		    contact1++)
		  for(contact2  = potential_waters[neighbour_indices[j]].contacts.begin();
		      contact2 != potential_waters[neighbour_indices[j]].contacts.end();
		      contact2++)
		    if(*contact1 == *contact2)
		      shared_contact = true;
	      
	      if(!shared_contact) //merge, if waters do not share a contact
		// if(fabs((*potential_water1)->water_atom.x  - potential_water2->water_atom.x)<water_join_dist)
		// if(fabs(potential_water1->water_atom.y  - potential_water2->water_atom.y)<water_join_dist)
		//    if(fabs(potential_water1->water_atom.z  - potential_water2->water_atom.z)<water_join_dist)
		if(dist.calculate(&potential_waters[indices[i]].water_atom,
				  &potential_waters[neighbour_indices[j]].water_atom, false) < water_join_dist_sq) //water atoms must be closer than join_dist
		  
		  {
		    
		    merge(&potential_waters[indices[i]],&potential_waters[neighbour_indices[j]]);
		  }
	      
	    }

  

  //cout<<"done"<<endl;
  /*
  for(potential_water1  = potential_waters.begin(); 
      potential_water1 != potential_waters.end(); 
      potential_water1++)
    if(potential_water1->active)
      for(potential_water2  = potential_water1+1; 
	  potential_water2 != potential_waters.end();
	  potential_water2++)
	if(potential_water2->active)
	  if(abs(potential_water1->water_atom.bx  - potential_water2->water_atom.bx)<2)
	    if(abs(potential_water1->water_atom.by  - potential_water2->water_atom.by)<2)
	      if(abs(potential_water1->water_atom.bz  - potential_water2->water_atom.bz)<2)
	    {
	   
	    shared_contact = false; //first check if waters have same primary contact
	    if(potential_water1->contacts[0] == potential_water2->contacts[0])
	      shared_contact = true;

	    if(!shared_contact) //if not go through all contacts
	      for(contact1  = potential_water1->contacts.begin();
		  contact1 != potential_water1->contacts.end();
		  contact1++)
		for(contact2  = potential_water2->contacts.begin();
		    contact2 != potential_water2->contacts.end();
		    contact2++)
		  if(*contact1 == *contact2)
		    shared_contact = true;
	    
	    if(!shared_contact) //merge, if waters do not share a contact
	      if(fabs(potential_water1->water_atom.x  - potential_water2->water_atom.x)<water_join_dist)
		 if(fabs(potential_water1->water_atom.y  - potential_water2->water_atom.y)<water_join_dist)
		    if(fabs(potential_water1->water_atom.z  - potential_water2->water_atom.z)<water_join_dist)
		       if(dist.calculate(&potential_water1->water_atom,
					 &potential_water2->water_atom, false) < water_join_dist_sq) //water atoms must be closer than join_dist
		
			 {
			   
			   merge(&*potential_water1,&*potential_water2);
			 }
		       
	  }

  */
}

void Water_Bridges::merge(potential_water_struct * potential_water1, potential_water_struct * potential_water2)
{
  // set joined water contact
  potential_water1->water_atom = Atom(Vector((potential_water1->water_atom.x + potential_water2->water_atom.x)/2.0,
					     (potential_water1->water_atom.y + potential_water2->water_atom.y)/2.0,
					     (potential_water1->water_atom.z + potential_water2->water_atom.z)/2.0),
				      potential_water1->water_atom.residue_no);
  
  // share contacts
  potential_water1->contacts.insert (potential_water1->contacts.end(), 
				     potential_water2->contacts.begin(),
				     potential_water2->contacts.end()); 
  
  // delete second contact
  potential_water2->active=false;
 
  
  return;

}


void Water_Bridges::remove_water_bumps(vector<int> indices, vector<int> neighbour_indices)
{
  for(unsigned int i=0; i<indices.size(); i++)
    if(potential_waters[indices[i]].active)
      for(unsigned int j=0; j<neighbour_indices.size(); j++)
	if(indices[i] != neighbour_indices[j])
	  if(potential_waters[neighbour_indices[j]].active)
	    if(dist.calculate(&potential_waters[indices[i]].water_atom, 
			      &potential_waters[neighbour_indices[j]].water_atom, false) < protein_bump_dist_sq )
	      {
		
		// delete the water with the fewest contacts
		if (potential_waters[indices[i]].contacts.size() > potential_waters[neighbour_indices[j]].contacts.size())
		  potential_waters[neighbour_indices[j]].active=false;
		else if (potential_waters[neighbour_indices[j]].contacts.size() > potential_waters[indices[i]].contacts.size())
		  potential_waters[indices[i]].active=false;
		else
		  potential_waters[indices[i]].active=false;
		

	      }
}

void Water_Bridges::set_parameters()
{
  same_residue_factor = 0.0;
  neighbour_residue_factor = 0.5;
  
  //ref: Guerois,Nielsen and Serrano: Predictig Changes in the stability of Proteins and Protein Complexes: 
  //                                  A study of more than 1000 mutations.
  water_bridge_energy = -1.4 * 4.184; //kJ/mol

  optimal_hb_dist = 2.85; //optimal hb dist  - Taylor & Kennard 1984
  min_hb_dist = 2.0;

  optimal_hb_angles.push_back(120);
  optimal_hb_angles.push_back(140);
  optimal_hb_angles.push_back(160);

  protein_bump_dist = 3.4;
  water_join_dist = 1.0;

  include_surf_factor = 1;

  //ref: Guerois,Nielsen and Serrano: Predictig Changes in the stability of Proteins and Protein Complexes: 
  //                                  A study of more than 1000 mutations.
  water_desolvation_penalty = 0.0;////////////////////////////////////////////////1.8*4.184;//kJ/mol 

  float step = 0.5;// Aangstroem

  for(unsigned int i=0; i<optimal_hb_angles.size(); i++)
    {
      d_along.push_back( optimal_hb_dist * fabs(cos(optimal_hb_angles[i] * 3.14159265/180)));
      d_perp.push_back( optimal_hb_dist * fabs(sin(optimal_hb_angles[i] * 3.14159265/180)));

      steps.push_back(step/(optimal_hb_dist*sin(3.14159265 - optimal_hb_angles[i] * 3.14159265/180))*180/3.14159265);
    }
  
  protein_bump_dist_sq = protein_bump_dist*protein_bump_dist;
  water_join_dist_sq = water_join_dist*water_join_dist;
  min_hb_dist_sq = min_hb_dist*min_hb_dist;
  optimal_hb_dist_sq = optimal_hb_dist*optimal_hb_dist;
  double_optimal_hb_dist_sq = 4*optimal_hb_dist*optimal_hb_dist;
}


bool Water_Bridges::no_heavy_atoms_bonded(Atom * atom)
{

  bool res = true;
  for(int i=0; i< atom->bonded_atoms.size(); i++)
    if (atom->bonded_atoms[i]->element != "H")
      res = false;


  if(res)
    {
      cout<<"Warning: No heavy bonded atoms for ";
      atom->display();

    }
  return res;
}

void Water_Bridges::writeDescription(FILE * reporter)
{
  /*** write out to out file **/
  version =string("$Id: water_bridges.cpp 18 2005-11-15 09:56:08Z chresten $");

  description = 
     string("");

  fprintf(reporter,"%s\n",version.c_str());	
  fprintf(reporter,"%s\n",description.c_str());
}
