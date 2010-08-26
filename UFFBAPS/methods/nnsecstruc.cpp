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
/*
$Id: nnsecstruc.cpp 6326 2010-08-26 16:08:21Z nielsen $
*/

#include "nnsecstruc.h"


Nnsecstruc::Nnsecstruc(FILE * reporter, vector<SoupObject *> trainSet, vector<SoupObject *> testSet):Method(reporter)
{
  /*** write out to out file **/
  writeDescription(reporter);
 
  setup();
 
  train(trainSet);

  test(testSet);

}


Nnsecstruc::Nnsecstruc(FILE * reporter, vector<SoupObject *> testSet):Method(reporter)
{
  /*** write out to out file **/
  writeDescription(reporter);
 
  setup();
  /*** use hard coded weights ***/
  setWeights();
  test(testSet);
}


Nnsecstruc::~Nnsecstruc(){}



void Nnsecstruc::train(vector<SoupObject *> trainSet)
{
  /*** read secondary structures from files ***/
  for(int testno=0; testno< trainSet.size(); testno++)
    {
      /*** casting soupobject to proteinchain to make residues available ***/
      ProteinChain * current = static_cast<ProteinChain *>(trainSet[testno]);
      
      /*** read secondary structure information form original file ***/
      string orifile = string(current->identify(), 0, 8);
      
      ifstream file(orifile.c_str());
      string line;
      helixResidues.clear();
      sheetResidues.clear();
	  
      if(!file.is_open())
	{
	  printf("Original file %s not found!\n",orifile.c_str());
	  exit(0);
	}
      
      while(!file.eof())
	{
	  getline(file, line);
	  string id = string(line, 0,6);
	  /*** reading info on helices ***/
	  if(id == "HELIX ")
	    {
	      int start = atoi(string(line,21,4).c_str());
	      int end   = atoi(string(line,33,4).c_str());
	      
	      int no = start;
	      while(no<=end)
		{
		  helixResidues.push_back(no);
		  no++;
		}
	    }
	  /*** reading info on sheets ***/
	  if(id == "SHEET " )
	    {
	      int start = atoi(string(line,22,4).c_str());
	      int end   = atoi(string(line,33,4).c_str());
	      
	      int no = start;
	      while(no<=end)
		{
		  sheetResidues.push_back(no);
		  no++;
		}	  
	    }
	}
      helixLibrary.push_back(helixResidues);
      sheetLibrary.push_back(sheetResidues);
    }



  /*** for each chain in training set ***/
  for(int run =0; run < 500; run++)
    {
      for(int testno=0; testno< trainSet.size(); testno++)
	{
	  /*** casting soupobject to proteinchain to make residues available ***/
	  ProteinChain * current = static_cast<ProteinChain *>(trainSet[testno]);
	  
	  helixResidues.clear();
	  sheetResidues.clear();

	  helixResidues = helixLibrary[testno];
	  sheetResidues = sheetLibrary[testno];

	  /*** start the training ***/
	  vector <double> temp;
	  temp.assign(2,0);
	  predictions.assign(current->residues.size(),temp);
	
	  for(int res=0; res< current->residues.size();res++)
	    {
	
	      vector<vector<double> > output = compute(current->residues, res);
	      vector<double> answer = getAnswer(res);

	      predictions[res][0] = output[0][1];
	      predictions[res][1] = output[1][1];

	      /*** errors on output units ***/
	  
	      errors[0][1] = output[0][1]*(1 - output[0][1])*(answer[0] - output[0][1]);
	      errors[1][1] = output[1][1]*(1 - output[1][1])*(answer[1] - output[1][1]);
	 
	      /*** errors on hidden units ***/
	      
	      errors[0][0] = 
		output[0][0]*(1 - output[0][0])* (weight12[0]*errors[0][1]+weight22[0]*errors[1][1]);
	      
	      errors[1][0] = 
		output[1][0]*(1 - output[1][0])* (weight12[1]*errors[0][1]+weight22[1]*errors[1][1]);
	      
	      
	      
	      /*** update weights to output units ***/
	      for(int i=0; i<2;i++)
		{
		  weight12[i] += learningRate * errors[0][1] * output[i][0];
		  weight22[i] += learningRate * errors[1][1] * output[i][0];
		}
	      

	      /*** update weights to hidden units ***/
	      for(int p=0; p<17;p++)
		for(int a=0; a<21;a++)
		  {
		    weight11[p][a] += learningRate * errors[0][0] * input[p][a];
		    weight21[p][a] += learningRate * errors[1][0] * input[p][a];
		  }
	      
	    } //loop over residues
	  printf("********* For %s\n",current->identify().c_str());
	  count(current->residues.size());	
	
	}//loop over trainingset


      flushHistogram();
      printf("After run %d I got %d right and %d wrong - %.3f \%\n", run, allright, allwrong, 100*((double) allright) /(allright + allwrong));
      
      allright =0;
      allwrong =0;
    } //loop over  500 runs
  

  printf("--------------------------------------------------------------\n");
  printf("Training done!\n");
  printf("--------------------------------------------------------------\n");
  
  for(int i=0;i<17;i++)
    for(int j=0;j<21;j++)
      printf("weight11[%d][%d] = %f;\n", i,j,weight11[i][j]);

   printf("--------------------------------------------------------------\n");

  for(int i=0;i<17;i++)
    for(int j=0;j<21;j++)
      printf("weight21[%d][%d] = %f;\n", i,j,weight21[i][j]);

  printf("--------------------------------------------------------------\n");
  for(int i=0;i<2;i++)
      printf("weight12[%d] = %f;\n",i,weight12[i]);

  printf("--------------------------------------------------------------\n");
  for(int i=0;i<2;i++)
      printf("weight22[%d] = %f;\n",i,weight22[i]);


}
void Nnsecstruc::test(vector<SoupObject *> testSet)
{


  allright += 0;
  allwrong += 0;

  for(int testno=0; testno< testSet.size(); testno++)
    {
      /*** casting soupobject to proteinchain to make residues available ***/
      ProteinChain * current = static_cast<ProteinChain *>(testSet[testno]);
	  
      /*** read secondary structure information form original file ***/
      string orifile = string(current->identify(), 0, 8);
	  
      ifstream file(orifile.c_str());
      string line;
      helixResidues.clear();
      sheetResidues.clear();
	  
      if(!file.is_open())
	{
	  printf("Original file %s not found!\n",orifile.c_str());
	  exit(0);
	}
      
      while(!file.eof())
	{
	  getline(file, line);
	  string id = string(line, 0,6);
	  /*** reading info on helices ***/
	  if(id == "HELIX ")
	    {
	      int start = atoi(string(line,21,4).c_str());
	      int end   = atoi(string(line,33,4).c_str());
	      
	      int no = start;
	      while(no<=end)
		{
		  helixResidues.push_back(no);
		  no++;
		}
	    }
	  /*** reading info on sheets ***/
	  if(id == "SHEET " )
	    {
	      int start = atoi(string(line,22,4).c_str());
	      int end   = atoi(string(line,33,4).c_str());
	      
	      int no = start;
	      while(no<=end)
		{
		  sheetResidues.push_back(no);
		  no++;
		}	  
	    }
	}
      
      /*** start the testing ***/
      vector <double> temp;
      temp.assign(2,0);
      predictions.assign(current->residues.size(),temp);
	
      for(int res=0; res< current->residues.size();res++)
	{
	  vector<vector<double> > output = compute(current->residues, res);
	  predictions[res][0] = output[0][1];
	  predictions[res][1] = output[1][1];
	 
	} //loop over residues
      
      /*** counting right predictions ***/
      printf("*********TEST*** For %s\n",current->identify().c_str());
      count(current->residues.size());
    }//loop over trainingset
  
      flushHistogram();
  printf("In testing I got %d right and %d wrong - %.3f \%\n", allright, allwrong,
	 100*((double) allright) /(allright + allwrong));
  
  
}

void Nnsecstruc::flushHistogram()
{
  int tot=0;
  for(int i = 0; i< histogram.size()-1;i++)
    {
      printf(" %4d  ",histogram[i]);
      tot +=histogram[i];
    } 
  printf("\n");
  for(int i = 0; i< histogram.size()-1;i++)
    printf("  %.2f ",((double) i)/(histogram.size()-1));
  printf("  sum: %d",tot);
  printf("\n");

  histogram.assign(11,0);

}

void Nnsecstruc::count(int size)
{
  /*** counting right predictions ***/
  int right =0;
  int wrong =0;

  for(int res=0; res< size;res++)
    {
      int no = (int) (predictions[res][0]*(histogram.size()-1));
      histogram[no]++;
      no = (int) (predictions[res][1]*(histogram.size()-1));
      histogram[no]++;
    }

  /*** setting output to 0 or 1 ***/
  double threshold = 0.37; 
  for(int res=0; res< size;res++)
    {
      vector<double> answer = getAnswer(res);
      
      /*** avoids 1 1 outputs ***/
      if(predictions[res][0] > threshold && predictions[res][1] > threshold)
	{
	  if(predictions[res][0] > predictions[res][1])
	    predictions[res][1] = 0;
	  else
	    predictions[res][0] = 0;
	}
	
      if(predictions[res][0] > threshold)
	predictions[res][0] = 1;
      else
	predictions[res][0] = 0;
      
      if(predictions[res][1] > threshold)
	predictions[res][1] = 1;
      else
	predictions[res][1] = 0;
    }
  
  /*** setting min helix lenght to 4 and min sheet lenght to 2 ***/

  int lenght =0;
  string current = "dummy";
  string now;
  for(int res=0; res< size;res++)
    {
      now = "coil";
      if(predictions[res][0] == 1  && predictions[res][1] == 0)
	now = "helix";
      if(predictions[res][0] == 0  && predictions[res][1] == 1)
	now = "sheet";

      /*** counting ***/
      if(now == "helix" && current == "helix")
	lenght++;

      if(now == "sheet" && current == "sheet")
	lenght++;

      /*** found end ***/
      if(current == "helix" && (now != "helix" || res == size-1))
	{
	  if(lenght < 4)
	    for(int pos = res-lenght; pos<res; pos++)
	      predictions[pos][0] = 0;

	  printf("Found helix of lenght %d and ending at %d\n",lenght,res);
	  lenght = 0;
	  current ="dummy";

	}
      if(current == "sheet" && (now != "sheet" || res == size-1))
	{
	  if(lenght < 2)
	    predictions[res-1][1] = 0;

	  printf("Found sheet of lenght %d and ending at %d\n",lenght,res);
	  lenght = 0;
	  current ="dummy";

	}

      /*** new helix ***/
      if(current == "dummy" && now == "helix")
	{
	  current = "helix";
	  lenght++;
	}
      /*** new sheet ***/
      if(current == "dummy" && now == "sheet")
	{
	  current = "sheet";
	  lenght++;
	}
    }
  for(int res=0; res< size;res++)
    {
      vector<double> answer = getAnswer(res);

      if(predictions[res][0] == answer[0] && predictions[res][1] == answer[1])
	right ++;
      else
	wrong ++;
     
      printf(" residue %3d was found to be ",res+1);
      if(predictions[res][0]==0 && predictions[res][1] ==0)
	printf("coil  ");
      if(predictions[res][0]==1 && predictions[res][1] ==0)
	printf("helix ");
      if(predictions[res][0]==0 && predictions[res][1] ==1)
	printf("sheet ");
      if(predictions[res][0]==1 && predictions[res][1] ==1)
	printf("ERROR ");

      printf("and should be ");

      if(answer[0]==0 && answer[1] ==0)
	printf("coil  ");
      if(answer[0]==1 && answer[1] ==0)
	printf("helix ");
      if(answer[0]==0 && answer[1] ==1)
	printf("sheet ");
      if(answer[0]==1 && answer[1] ==1)
	printf("ERROR ");

      
      if(predictions[res][0] == answer[0] && predictions[res][1] == answer[1])
	printf("     TRUE\n");
      else
	printf("   FALSE\n");

    }
  
  allright += right;
  allwrong += wrong;
	
}


vector<vector<double> > Nnsecstruc::compute(vector<Residue> residues, int res)
{
  int min = 0;
  int max = residues.size()-1;

  vector<vector<double> > result;
  vector<double> temp;
  temp.assign(2,0);
  result.assign(2,temp);


  /*** assigning input ***/
  for(int p=0; p<17;p++)
    {
      int curA=0;
      int pos = res - 9 + p;
      if(pos >= min && pos <= max)
	  curA = conversion[residues[pos].residueType];

      for(int a=0; a<21;a++)
	{
	  if(a == curA)
	    input[p][a] = 1;
	  else
	    input[p][a] = 0;
	}
    }

  hiddenOutput.assign(2,0);
  output.assign(2,0);

  /*** hidden units ***/
  for(int p=0; p<17;p++)
    for(int a=0; a<21;a++)
      {
	hiddenOutput[0] +=  weight11[p][a]*input[p][a];
	hiddenOutput[1] +=  weight21[p][a]*input[p][a];
      }
 
  hiddenOutput[0] = sgm(hiddenOutput[0]);
  hiddenOutput[1] = sgm(hiddenOutput[1]);
    
  /*** output units ***/
  for(int i=0; i<2;i++)
    {
      output[0] += weight12[i]*hiddenOutput[i];
      output[1] += weight22[i]*hiddenOutput[i];
    }
  
  output[0] = sgm(output[0]);
  output[1] = sgm(output[1]);

  result[0][0] = hiddenOutput[0];
  result[1][0] = hiddenOutput[1];
  result[0][1] = output[0];
  result[1][1] = output[1];

  return result;
}

double Nnsecstruc::sgm(double in){ 

  double res = 1/(1+exp(-in));
  return res;
}



vector<double> Nnsecstruc::getAnswer(int res)
{
  /*** looks up right answer ***/
  vector<double> answer;
  answer.push_back(0.0);
  answer.push_back(0.0);

  vector<int>::iterator helix;
  vector<int>::iterator sheet;
  helix = find( helixResidues.begin(), helixResidues.end(), res+1);
  sheet = find( sheetResidues.begin(), sheetResidues.end(), res+1);
  
  if(helix != helixResidues.end()) 
    {
      answer[0] = 1;
      answer[1] = 0;
    }
  else if(sheet != sheetResidues.end()) 
    {
      answer[0] = 0;
      answer[1] = 1;
    }
  return answer;

}


void Nnsecstruc::setWeights()
{

  /*** hard coded weights from training on
       
  1hel.pdb 1lyz.pdb 1bin.pdb 1fox.pdb 1c2o.pdb 2rmb.pdb 1pin.pdb 1sbt.pdb 8paz.pdb 3cna.pdb
  
  After 500 epochs secondary structures of the training set 
  was predicted to an accurancy of 81.361 %
  and 31.735 % for the test set consisting of 
  
       1mbn.pdb 1eco.pdb 2lhb.pdb
 ***/

  //--------------------------------------------------------------
weight11[0][0] = -3.090259;
weight11[0][1] = 0.582326;
weight11[0][2] = 3.384177;
weight11[0][3] = 1.112229;
weight11[0][4] = -1.329713;
weight11[0][5] = 0.148475;
weight11[0][6] = -0.135143;
weight11[0][7] = 0.558443;
weight11[0][8] = 0.045844;
weight11[0][9] = 0.340577;
weight11[0][10] = -0.598584;
weight11[0][11] = 0.923774;
weight11[0][12] = -1.298490;
weight11[0][13] = 3.077046;
weight11[0][14] = 1.404438;
weight11[0][15] = 0.472689;
weight11[0][16] = 1.226007;
weight11[0][17] = -1.509851;
weight11[0][18] = 1.650671;
weight11[0][19] = -1.317717;
weight11[0][20] = 2.134825;
weight11[1][0] = 0.994087;
weight11[1][1] = 0.211237;
weight11[1][2] = 3.195187;
weight11[1][3] = -1.103046;
weight11[1][4] = 0.224610;
weight11[1][5] = 2.925375;
weight11[1][6] = -0.364939;
weight11[1][7] = 1.690317;
weight11[1][8] = 0.534843;
weight11[1][9] = -2.710820;
weight11[1][10] = -0.818636;
weight11[1][11] = -0.889716;
weight11[1][12] = 0.944254;
weight11[1][13] = 5.037066;
weight11[1][14] = -0.671358;
weight11[1][15] = -1.865603;
weight11[1][16] = 0.690117;
weight11[1][17] = -1.725079;
weight11[1][18] = 0.886722;
weight11[1][19] = 0.431505;
weight11[1][20] = 0.136039;
weight11[2][0] = 1.690355;
weight11[2][1] = 0.704626;
weight11[2][2] = 2.457044;
weight11[2][3] = -1.031129;
weight11[2][4] = 0.601045;
weight11[2][5] = 1.472524;
weight11[2][6] = 1.603196;
weight11[2][7] = 0.260407;
weight11[2][8] = 0.687439;
weight11[2][9] = -0.786652;
weight11[2][10] = -0.513427;
weight11[2][11] = -0.718206;
weight11[2][12] = -1.417368;
weight11[2][13] = 2.185347;
weight11[2][14] = 0.972766;
weight11[2][15] = 1.886466;
weight11[2][16] = 1.233092;
weight11[2][17] = 0.135125;
weight11[2][18] = -3.421985;
weight11[2][19] = -0.686401;
weight11[2][20] = 0.555602;
weight11[3][0] = 2.346079;
weight11[3][1] = 2.724004;
weight11[3][2] = -0.634658;
weight11[3][3] = -2.049344;
weight11[3][4] = 0.350631;
weight11[3][5] = 0.702732;
weight11[3][6] = -0.487913;
weight11[3][7] = 4.510212;
weight11[3][8] = 0.227727;
weight11[3][9] = -2.409876;
weight11[3][10] = 1.150445;
weight11[3][11] = 0.052304;
weight11[3][12] = -0.955627;
weight11[3][13] = 3.417404;
weight11[3][14] = 0.499547;
weight11[3][15] = -0.178469;
weight11[3][16] = 0.169701;
weight11[3][17] = -2.084778;
weight11[3][18] = 0.323055;
weight11[3][19] = 0.510877;
weight11[3][20] = -0.544531;
weight11[4][0] = 1.801418;
weight11[4][1] = 0.519947;
weight11[4][2] = -0.482402;
weight11[4][3] = 0.532297;
weight11[4][4] = -1.443826;
weight11[4][5] = 6.515918;
weight11[4][6] = 1.201611;
weight11[4][7] = 3.264743;
weight11[4][8] = -1.399210;
weight11[4][9] = -2.070385;
weight11[4][10] = 1.501882;
weight11[4][11] = -1.202155;
weight11[4][12] = -1.043879;
weight11[4][13] = -0.802226;
weight11[4][14] = 3.591280;
weight11[4][15] = 2.085022;
weight11[4][16] = 0.107650;
weight11[4][17] = -2.671504;
weight11[4][18] = 0.192337;
weight11[4][19] = 2.185699;
weight11[4][20] = -4.358623;
weight11[5][0] = 1.311913;
weight11[5][1] = -1.731577;
weight11[5][2] = -0.973313;
weight11[5][3] = 1.426780;
weight11[5][4] = -0.224562;
weight11[5][5] = 3.734266;
weight11[5][6] = 3.213252;
weight11[5][7] = 0.490666;
weight11[5][8] = 0.549693;
weight11[5][9] = -0.371382;
weight11[5][10] = -1.306612;
weight11[5][11] = 2.016388;
weight11[5][12] = 0.564701;
weight11[5][13] = -2.667261;
weight11[5][14] = -2.339995;
weight11[5][15] = 3.712840;
weight11[5][16] = 2.445449;
weight11[5][17] = -0.260549;
weight11[5][18] = -1.915181;
weight11[5][19] = 3.634946;
weight11[5][20] = -3.446106;
weight11[6][0] = 3.426413;
weight11[6][1] = 0.147669;
weight11[6][2] = -1.307138;
weight11[6][3] = 2.027536;
weight11[6][4] = -0.272817;
weight11[6][5] = 4.997473;
weight11[6][6] = -1.129589;
weight11[6][7] = 6.847147;
weight11[6][8] = -0.462658;
weight11[6][9] = -3.781806;
weight11[6][10] = -0.494443;
weight11[6][11] = 2.294823;
weight11[6][12] = -0.492538;
weight11[6][13] = -1.333544;
weight11[6][14] = 0.011640;
weight11[6][15] = 0.954572;
weight11[6][16] = 1.785742;
weight11[6][17] = -3.533745;
weight11[6][18] = 2.081557;
weight11[6][19] = 1.434668;
weight11[6][20] = -5.225764;
weight11[7][0] = 1.870564;
weight11[7][1] = -0.072330;
weight11[7][2] = 2.352753;
weight11[7][3] = 6.353479;
weight11[7][4] = 0.871845;
weight11[7][5] = 0.559156;
weight11[7][6] = -1.895241;
weight11[7][7] = 3.939189;
weight11[7][8] = 1.648977;
weight11[7][9] = -1.478340;
weight11[7][10] = -1.376450;
weight11[7][11] = -1.312244;
weight11[7][12] = 2.177676;
weight11[7][13] = -2.642296;
weight11[7][14] = -4.398132;
weight11[7][15] = 4.742258;
weight11[7][16] = 1.506922;
weight11[7][17] = -1.241245;
weight11[7][18] = -0.375727;
weight11[7][19] = 1.330844;
weight11[7][20] = -4.601244;
weight11[8][0] = 1.743841;
weight11[8][1] = 2.352011;
weight11[8][2] = -0.083456;
weight11[8][3] = 2.429772;
weight11[8][4] = 0.284140;
weight11[8][5] = 0.686157;
weight11[8][6] = 2.993964;
weight11[8][7] = 5.998774;
weight11[8][8] = 2.239329;
weight11[8][9] = 1.400644;
weight11[8][10] = -4.445012;
weight11[8][11] = -2.012751;
weight11[8][12] = 3.072553;
weight11[8][13] = -4.304324;
weight11[8][14] = -2.609563;
weight11[8][15] = 3.953325;
weight11[8][16] = 0.900116;
weight11[8][17] = -3.263164;
weight11[8][18] = 1.014929;
weight11[8][19] = -0.991391;
weight11[8][20] = -3.476941;
weight11[9][0] = 0.036553;
weight11[9][1] = 0.409406;
weight11[9][2] = -3.011980;
weight11[9][3] = 2.483105;
weight11[9][4] = 4.188640;
weight11[9][5] = 0.764521;
weight11[9][6] = 2.526197;
weight11[9][7] = 3.293499;
weight11[9][8] = -0.684294;
weight11[9][9] = 1.739073;
weight11[9][10] = -3.882041;
weight11[9][11] = -0.024423;
weight11[9][12] = 0.692123;
weight11[9][13] = 0.920893;
weight11[9][14] = -0.924049;
weight11[9][15] = 5.640838;
weight11[9][16] = -0.583995;
weight11[9][17] = -3.561612;
weight11[9][18] = 1.252178;
weight11[9][19] = -1.207240;
weight11[9][20] = -2.317211;
weight11[10][0] = 0.305150;
weight11[10][1] = 0.393536;
weight11[10][2] = -1.208951;
weight11[10][3] = 3.593425;
weight11[10][4] = 3.311005;
weight11[10][5] = 0.133322;
weight11[10][6] = 1.465858;
weight11[10][7] = 1.738484;
weight11[10][8] = 0.425899;
weight11[10][9] = 0.185173;
weight11[10][10] = -3.809399;
weight11[10][11] = 1.677942;
weight11[10][12] = 1.720219;
weight11[10][13] = 0.905902;
weight11[10][14] = -0.411360;
weight11[10][15] = 1.738273;
weight11[10][16] = 0.132145;
weight11[10][17] = -3.972974;
weight11[10][18] = -0.431957;
weight11[10][19] = 0.813453;
weight11[10][20] = -1.049806;
weight11[11][0] = 0.001176;
weight11[11][1] = -2.701609;
weight11[11][2] = -3.342375;
weight11[11][3] = 4.813019;
weight11[11][4] = 2.600371;
weight11[11][5] = -0.595566;
weight11[11][6] = 4.184691;
weight11[11][7] = 3.004326;
weight11[11][8] = 1.139309;
weight11[11][9] = 2.868655;
weight11[11][10] = -1.691496;
weight11[11][11] = 0.110705;
weight11[11][12] = 1.202214;
weight11[11][13] = -1.612476;
weight11[11][14] = -2.748758;
weight11[11][15] = 4.125676;
weight11[11][16] = 1.359071;
weight11[11][17] = -1.293897;
weight11[11][18] = 0.086506;
weight11[11][19] = -4.354483;
weight11[11][20] = 0.655337;
weight11[12][0] = -2.588354;
weight11[12][1] = -1.229612;
weight11[12][2] = -2.994931;
weight11[12][3] = 0.522392;
weight11[12][4] = 1.999718;
weight11[12][5] = -0.659925;
weight11[12][6] = 2.388273;
weight11[12][7] = 5.321957;
weight11[12][8] = 1.753073;
weight11[12][9] = -0.371407;
weight11[12][10] = -2.884066;
weight11[12][11] = 0.750212;
weight11[12][12] = 2.955533;
weight11[12][13] = -0.156960;
weight11[12][14] = -0.347836;
weight11[12][15] = 2.551083;
weight11[12][16] = -0.048803;
weight11[12][17] = -1.625359;
weight11[12][18] = 0.309434;
weight11[12][19] = 0.067811;
weight11[12][20] = 2.390474;
weight11[13][0] = 0.999523;
weight11[13][1] = 0.239384;
weight11[13][2] = 1.802930;
weight11[13][3] = 0.123948;
weight11[13][4] = -2.553351;
weight11[13][5] = -3.274671;
weight11[13][6] = 1.186268;
weight11[13][7] = 2.102769;
weight11[13][8] = 1.468599;
weight11[13][9] = 1.835751;
weight11[13][10] = -1.878601;
weight11[13][11] = 3.692495;
weight11[13][12] = 1.022887;
weight11[13][13] = -0.454726;
weight11[13][14] = -2.228016;
weight11[13][15] = -0.343544;
weight11[13][16] = 0.807343;
weight11[13][17] = 1.027796;
weight11[13][18] = 1.311570;
weight11[13][19] = 0.547203;
weight11[13][20] = 0.379244;
weight11[14][0] = 0.143282;
weight11[14][1] = 0.770079;
weight11[14][2] = 0.000912;
weight11[14][3] = 0.905757;
weight11[14][4] = -1.972305;
weight11[14][5] = 1.852059;
weight11[14][6] = 2.724581;
weight11[14][7] = 2.283769;
weight11[14][8] = 1.976362;
weight11[14][9] = 2.442993;
weight11[14][10] = 0.985910;
weight11[14][11] = 4.144940;
weight11[14][12] = 1.887845;
weight11[14][13] = 1.262217;
weight11[14][14] = -3.007869;
weight11[14][15] = -0.566915;
weight11[14][16] = -0.654267;
weight11[14][17] = -2.361677;
weight11[14][18] = -3.891897;
weight11[14][19] = -0.649766;
weight11[14][20] = -0.442354;
weight11[15][0] = 1.668025;
weight11[15][1] = -0.257514;
weight11[15][2] = 4.101538;
weight11[15][3] = 0.679332;
weight11[15][4] = -0.196485;
weight11[15][5] = 1.336531;
weight11[15][6] = 0.139540;
weight11[15][7] = 0.891322;
weight11[15][8] = 1.678966;
weight11[15][9] = -1.146331;
weight11[15][10] = -0.106005;
weight11[15][11] = 0.829624;
weight11[15][12] = -0.357066;
weight11[15][13] = 2.788393;
weight11[15][14] = -0.164137;
weight11[15][15] = 0.371077;
weight11[15][16] = -1.142180;
weight11[15][17] = -5.427409;
weight11[15][18] = 1.324031;
weight11[15][19] = -0.484290;
weight11[15][20] = 1.424873;
weight11[16][0] = 6.634439;
weight11[16][1] = 0.060673;
weight11[16][2] = 3.061142;
weight11[16][3] = -2.869014;
weight11[16][4] = 1.074767;
weight11[16][5] = -0.130832;
weight11[16][6] = -0.799358;
weight11[16][7] = 2.091027;
weight11[16][8] = 1.421561;
weight11[16][9] = -1.701648;
weight11[16][10] = 0.088040;
weight11[16][11] = 1.581807;
weight11[16][12] = 1.459983;
weight11[16][13] = 0.944045;
weight11[16][14] = 1.677603;
weight11[16][15] = -0.113629;
weight11[16][16] = -2.794628;
weight11[16][17] = -5.284545;
weight11[16][18] = 1.484757;
weight11[16][19] = -2.589463;
weight11[16][20] = 2.508590;
//--------------------------------------------------------------
weight21[0][0] = 2.817606;
weight21[0][1] = -2.052392;
weight21[0][2] = -1.281897;
weight21[0][3] = 3.299582;
weight21[0][4] = -0.159061;
weight21[0][5] = -0.331413;
weight21[0][6] = -1.976677;
weight21[0][7] = -0.271722;
weight21[0][8] = -0.094835;
weight21[0][9] = 6.586731;
weight21[0][10] = -0.748769;
weight21[0][11] = -1.120206;
weight21[0][12] = -1.063683;
weight21[0][13] = 1.035779;
weight21[0][14] = -1.245905;
weight21[0][15] = -1.186427;
weight21[0][16] = 0.710176;
weight21[0][17] = 0.877962;
weight21[0][18] = -0.793901;
weight21[0][19] = 0.753862;
weight21[0][20] = 1.935573;
weight21[1][0] = 0.040682;
weight21[1][1] = -2.632165;
weight21[1][2] = 0.401545;
weight21[1][3] = 1.826642;
weight21[1][4] = -0.092726;
weight21[1][5] = -0.100610;
weight21[1][6] = -4.068076;
weight21[1][7] = 2.253051;
weight21[1][8] = 2.046768;
weight21[1][9] = 0.542945;
weight21[1][10] = -1.480154;
weight21[1][11] = -2.102095;
weight21[1][12] = 4.510824;
weight21[1][13] = 3.106664;
weight21[1][14] = -4.878132;
weight21[1][15] = 1.571236;
weight21[1][16] = 2.413578;
weight21[1][17] = 1.728553;
weight21[1][18] = -0.200440;
weight21[1][19] = -0.766628;
weight21[1][20] = 1.555897;
weight21[2][0] = 0.544208;
weight21[2][1] = -3.152961;
weight21[2][2] = 0.307246;
weight21[2][3] = -1.177369;
weight21[2][4] = -0.156026;
weight21[2][5] = 0.424477;
weight21[2][6] = -0.363876;
weight21[2][7] = -1.816785;
weight21[2][8] = 1.235051;
weight21[2][9] = 0.142252;
weight21[2][10] = 1.857899;
weight21[2][11] = -1.251830;
weight21[2][12] = 0.828041;
weight21[2][13] = 3.419271;
weight21[2][14] = -2.215713;
weight21[2][15] = 3.877330;
weight21[2][16] = 1.504681;
weight21[2][17] = 0.771164;
weight21[2][18] = 0.325983;
weight21[2][19] = -0.054218;
weight21[2][20] = 0.706381;
weight21[3][0] = -0.461103;
weight21[3][1] = -2.223877;
weight21[3][2] = -1.364471;
weight21[3][3] = -1.479203;
weight21[3][4] = 0.186595;
weight21[3][5] = 0.631642;
weight21[3][6] = -1.106185;
weight21[3][7] = 0.672879;
weight21[3][8] = 1.844746;
weight21[3][9] = 3.013193;
weight21[3][10] = -0.408601;
weight21[3][11] = 2.600206;
weight21[3][12] = 0.727291;
weight21[3][13] = 2.433624;
weight21[3][14] = -2.271269;
weight21[3][15] = 1.563888;
weight21[3][16] = 0.717115;
weight21[3][17] = -0.125875;
weight21[3][18] = -1.373733;
weight21[3][19] = 0.361266;
weight21[3][20] = 1.563010;
weight21[4][0] = 0.423538;
weight21[4][1] = -0.644712;
weight21[4][2] = -0.401970;
weight21[4][3] = -0.701553;
weight21[4][4] = -0.500759;
weight21[4][5] = -1.526742;
weight21[4][6] = -0.007830;
weight21[4][7] = 0.688407;
weight21[4][8] = 0.806056;
weight21[4][9] = 5.778394;
weight21[4][10] = 3.436401;
weight21[4][11] = 0.709755;
weight21[4][12] = 0.922227;
weight21[4][13] = -1.348756;
weight21[4][14] = 0.206351;
weight21[4][15] = 3.195344;
weight21[4][16] = 0.049532;
weight21[4][17] = -1.655100;
weight21[4][18] = -2.962582;
weight21[4][19] = -2.152274;
weight21[4][20] = 1.416055;
weight21[5][0] = 5.459678;
weight21[5][1] = -3.064859;
weight21[5][2] = 0.152673;
weight21[5][3] = 1.616022;
weight21[5][4] = 0.924462;
weight21[5][5] = -1.543359;
weight21[5][6] = -2.057853;
weight21[5][7] = 0.480862;
weight21[5][8] = 2.344273;
weight21[5][9] = 1.284859;
weight21[5][10] = 1.399339;
weight21[5][11] = 1.651500;
weight21[5][12] = 1.409614;
weight21[5][13] = -3.320648;
weight21[5][14] = 0.069190;
weight21[5][15] = 3.236667;
weight21[5][16] = 2.556415;
weight21[5][17] = -2.914857;
weight21[5][18] = -2.514312;
weight21[5][19] = -1.044449;
weight21[5][20] = -0.737997;
weight21[6][0] = 3.501749;
weight21[6][1] = -3.385958;
weight21[6][2] = 0.655955;
weight21[6][3] = 2.349225;
weight21[6][4] = -0.056060;
weight21[6][5] = -0.594361;
weight21[6][6] = 1.422012;
weight21[6][7] = 0.895071;
weight21[6][8] = 2.750307;
weight21[6][9] = -1.508053;
weight21[6][10] = 0.159557;
weight21[6][11] = -0.874631;
weight21[6][12] = -0.550440;
weight21[6][13] = -3.815421;
weight21[6][14] = 0.597501;
weight21[6][15] = 3.426240;
weight21[6][16] = 3.618259;
weight21[6][17] = 0.384914;
weight21[6][18] = -0.785805;
weight21[6][19] = -1.637158;
weight21[6][20] = -1.175454;
weight21[7][0] = 0.473159;
weight21[7][1] = -1.722970;
weight21[7][2] = -0.203067;
weight21[7][3] = 2.467033;
weight21[7][4] = 1.553283;
weight21[7][5] = -2.156730;
weight21[7][6] = -0.944161;
weight21[7][7] = 2.558644;
weight21[7][8] = 2.325234;
weight21[7][9] = -2.564389;
weight21[7][10] = -0.015767;
weight21[7][11] = -3.288995;
weight21[7][12] = 3.143304;
weight21[7][13] = -2.754729;
weight21[7][14] = -2.358201;
weight21[7][15] = 3.638265;
weight21[7][16] = 4.126122;
weight21[7][17] = 2.802913;
weight21[7][18] = 0.840066;
weight21[7][19] = -0.434080;
weight21[7][20] = -1.860943;
weight21[8][0] = 0.185681;
weight21[8][1] = -1.391704;
weight21[8][2] = -1.075869;
weight21[8][3] = 2.516323;
weight21[8][4] = 1.395332;
weight21[8][5] = -0.881011;
weight21[8][6] = 1.226997;
weight21[8][7] = 1.745713;
weight21[8][8] = 5.422809;
weight21[8][9] = -2.898758;
weight21[8][10] = -1.129642;
weight21[8][11] = -4.403765;
weight21[8][12] = 2.204747;
weight21[8][13] = -2.051233;
weight21[8][14] = 1.203224;
weight21[8][15] = 3.010523;
weight21[8][16] = 4.546530;
weight21[8][17] = 1.967901;
weight21[8][18] = -2.803356;
weight21[8][19] = -2.691328;
weight21[8][20] = -0.800772;
weight21[9][0] = 0.022101;
weight21[9][1] = -3.688034;
weight21[9][2] = -2.422383;
weight21[9][3] = 2.030883;
weight21[9][4] = 2.513373;
weight21[9][5] = -3.243358;
weight21[9][6] = 0.141632;
weight21[9][7] = 2.984791;
weight21[9][8] = 1.800255;
weight21[9][9] = -3.341119;
weight21[9][10] = -1.375326;
weight21[9][11] = -2.852268;
weight21[9][12] = 2.841818;
weight21[9][13] = 2.346820;
weight21[9][14] = 2.352169;
weight21[9][15] = 3.468261;
weight21[9][16] = 1.328785;
weight21[9][17] = 2.027128;
weight21[9][18] = -1.622118;
weight21[9][19] = 1.073752;
weight21[9][20] = -0.935697;
weight21[10][0] = 0.622177;
weight21[10][1] = -5.615015;
weight21[10][2] = -1.356142;
weight21[10][3] = 3.971934;
weight21[10][4] = 2.030328;
weight21[10][5] = 0.147393;
weight21[10][6] = -1.788087;
weight21[10][7] = -1.279866;
weight21[10][8] = -0.490840;
weight21[10][9] = -3.094372;
weight21[10][10] = -0.645763;
weight21[10][11] = -1.143863;
weight21[10][12] = 1.057223;
weight21[10][13] = 3.989556;
weight21[10][14] = 0.795007;
weight21[10][15] = 2.212914;
weight21[10][16] = 0.451394;
weight21[10][17] = 4.685465;
weight21[10][18] = -2.062918;
weight21[10][19] = 4.544182;
weight21[10][20] = -1.447049;
weight21[11][0] = 3.454818;
weight21[11][1] = -5.337424;
weight21[11][2] = -3.368021;
weight21[11][3] = 4.264655;
weight21[11][4] = 4.068051;
weight21[11][5] = 0.072205;
weight21[11][6] = -5.087916;
weight21[11][7] = 0.188669;
weight21[11][8] = 1.916310;
weight21[11][9] = -5.445365;
weight21[11][10] = 1.109763;
weight21[11][11] = -0.516861;
weight21[11][12] = 0.086177;
weight21[11][13] = -0.831612;
weight21[11][14] = 3.763379;
weight21[11][15] = 2.007552;
weight21[11][16] = 0.224868;
weight21[11][17] = 4.899371;
weight21[11][18] = -0.106516;
weight21[11][19] = 1.353150;
weight21[11][20] = -1.257028;
weight21[12][0] = -1.334164;
weight21[12][1] = -3.719203;
weight21[12][2] = -2.803681;
weight21[12][3] = 2.837736;
weight21[12][4] = 3.041131;
weight21[12][5] = -2.111837;
weight21[12][6] = -0.674536;
weight21[12][7] = -1.005260;
weight21[12][8] = 2.571605;
weight21[12][9] = -4.393970;
weight21[12][10] = 1.161528;
weight21[12][11] = -0.714420;
weight21[12][12] = -1.385390;
weight21[12][13] = 0.159420;
weight21[12][14] = 4.790769;
weight21[12][15] = 2.985951;
weight21[12][16] = -0.213736;
weight21[12][17] = 1.573996;
weight21[12][18] = 2.198680;
weight21[12][19] = 1.946722;
weight21[12][20] = 0.894217;
weight21[13][0] = -0.790578;
weight21[13][1] = -3.240383;
weight21[13][2] = 0.612064;
weight21[13][3] = -1.219170;
weight21[13][4] = -0.414801;
weight21[13][5] = -0.074606;
weight21[13][6] = -0.321368;
weight21[13][7] = 0.122255;
weight21[13][8] = 2.642186;
weight21[13][9] = 3.516234;
weight21[13][10] = 0.110241;
weight21[13][11] = 0.602762;
weight21[13][12] = -2.076945;
weight21[13][13] = -3.347818;
weight21[13][14] = 0.640151;
weight21[13][15] = 2.235488;
weight21[13][16] = -1.088033;
weight21[13][17] = 4.142097;
weight21[13][18] = 3.220500;
weight21[13][19] = 3.810032;
weight21[13][20] = -3.436157;
weight21[14][0] = 3.298641;
weight21[14][1] = -2.416247;
weight21[14][2] = 1.294433;
weight21[14][3] = -2.917559;
weight21[14][4] = 0.129899;
weight21[14][5] = -3.202417;
weight21[14][6] = 1.368755;
weight21[14][7] = -1.380207;
weight21[14][8] = 3.504350;
weight21[14][9] = 1.900358;
weight21[14][10] = 1.821769;
weight21[14][11] = 2.595568;
weight21[14][12] = -0.348518;
weight21[14][13] = 1.070012;
weight21[14][14] = -2.801255;
weight21[14][15] = -1.328096;
weight21[14][16] = -0.080831;
weight21[14][17] = -1.051370;
weight21[14][18] = 1.762803;
weight21[14][19] = 1.412746;
weight21[14][20] = 0.948808;
weight21[15][0] = 1.916838;
weight21[15][1] = -1.515888;
weight21[15][2] = 6.375990;
weight21[15][3] = -2.411804;
weight21[15][4] = -1.624136;
weight21[15][5] = 2.746416;
weight21[15][6] = -1.563346;
weight21[15][7] = -2.636768;
weight21[15][8] = 4.362746;
weight21[15][9] = -0.725536;
weight21[15][10] = -2.374808;
weight21[15][11] = -0.298631;
weight21[15][12] = -3.186807;
weight21[15][13] = 2.449993;
weight21[15][14] = -1.811456;
weight21[15][15] = 3.703146;
weight21[15][16] = 1.945065;
weight21[15][17] = -4.551173;
weight21[15][18] = 0.971333;
weight21[15][19] = 3.794754;
weight21[15][20] = -0.169397;
weight21[16][0] = 2.617556;
weight21[16][1] = -1.152341;
weight21[16][2] = 3.819802;
weight21[16][3] = -4.728518;
weight21[16][4] = -0.173191;
weight21[16][5] = -1.359004;
weight21[16][6] = -0.389692;
weight21[16][7] = -2.384848;
weight21[16][8] = 2.570403;
weight21[16][9] = 0.026527;
weight21[16][10] = 1.881307;
weight21[16][11] = 1.845768;
weight21[16][12] = 0.133888;
weight21[16][13] = 1.694013;
weight21[16][14] = -1.349459;
weight21[16][15] = -0.534807;
weight21[16][16] = -0.600260;
weight21[16][17] = -1.800354;
weight21[16][18] = 2.193143;
weight21[16][19] = 1.061845;
weight21[16][20] = 2.011806;
//--------------------------------------------------------------
weight12[0] = 14.489693;
weight12[1] = -17.222624;
//--------------------------------------------------------------
weight22[0] = -24.192253;
weight22[1] = 11.579582;


}

void Nnsecstruc::setup()
{
  allright =0;
  allwrong=0;
  int frameSize = 17;
  int noHidden = 2;
  learningRate = 0.15;
  histogram.assign(11,0);

  hiddenOutput.assign(2,0);
  output.assign(2,0);

  vector<int> temp;
  temp.assign(21,0);
  input.assign(frameSize,temp);
  


  vector<double> temp2;
  temp2.assign(21,0);
  weight11.assign(frameSize,temp2);
  weight21.assign(frameSize,temp2);

  temp2.clear();
  temp2.assign(noHidden,0);
  errors.assign(noHidden,temp2);


  weight12.assign(noHidden,0);
  weight22.assign(noHidden,0);

  // srand(time(NULL));

  
  for(int i=0;i<frameSize;i++)
    for(int j=0;j<21;j++)
      {
	weight11[i][j] = randomNo();
	weight21[i][j] = randomNo();
      }

  for(int i=0;i<noHidden;i++)
    {
      weight12[i] =randomNo();
      weight22[i] =randomNo();
    }

  conversion[string("")]    = 0;
  conversion[string("ALA")] = 1;
  conversion[string("ARG")] = 2;
  conversion[string("ASN")] = 3;
  conversion[string("ASP")] = 4;
  conversion[string("CYS")] = 5;
  conversion[string("GLN")] = 6;
  conversion[string("GLU")] = 7;
  conversion[string("GLY")] = 8;
  conversion[string("HIS")] = 9;
  conversion[string("ILE")] = 10;
  conversion[string("LEU")] = 11;
  conversion[string("LYS")] = 12;
  conversion[string("MET")] = 13;
  conversion[string("PHE")] = 14;
  conversion[string("PRO")] = 15;
  conversion[string("SER")] = 16;
  conversion[string("THR")] = 17;
  conversion[string("TRP")] = 18;
  conversion[string("TYR")] = 19;
  conversion[string("VAL")] = 20;

}


double Nnsecstruc::randomNo()
{
  double min = -0.05;
  double max =  0.05;

  int randNo = rand();

  double norm = ((double)randNo)/((double)RAND_MAX);

  double res = min + norm*(max - min);
  
  if(res >max)
    printf( "Random number: min %f, max %f, res %f \n", min, max, res );	
 
  return res;
}



void Nnsecstruc::writeDescription(FILE * reporter)
{
  /*** write out to out file **/
  version =string("$Id: nnsecstruc.cpp 6326 2010-08-26 16:08:21Z nielsen $");

  description = 
     string("Neural network for prediction secondary structure elements.");

  fprintf(reporter,"%s\n",version.c_str());	
  fprintf(reporter,"%s\n",description.c_str());
}
