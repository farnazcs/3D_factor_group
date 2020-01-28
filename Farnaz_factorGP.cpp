
#include <iostream>
#include <math.h>
#include <vector>
#include "./eigen-git-mirror/Eigen/Dense"
#include "./eigen-git-mirror/Eigen/LU"
#include <fstream>



using namespace std;

//This function gets two matrixed with NumRows rows and NumCols columns and returns true if they are close enough to each other
bool Matrixes_are_equal(Eigen::MatrixXd M1,Eigen::MatrixXd M2, double delta=0.001, int NumRows=3,int NumCols=3)
{
	double diff=0;
	// we subtract two mutrixes element-wise and add the differences together.
	for (int i=0;i<NumRows;i++)
		for(int j=0;j<NumCols;j++)
			diff=diff+fabs(M1(i,j)-M2(i,j));
	//if difference is small, they are similar or equal
	if (diff<delta)
		return true;
	else
		return false;
}
// this function generates the unity matrix of size LattDim.
Eigen::MatrixXd Unity_matrix(int LattDim=3)
{
	Eigen::MatrixXd UnityMat(LattDim,LattDim);
	//Setting elements (i,i) to 1 and the rest to zero.
	for (int i=0;i<LattDim;i++)
		for(int j=0;j<LattDim;j++)
			if(i==j)
				UnityMat(i,j)=1;
			else
				UnityMat(i,j)=0;
	return UnityMat;
}

//This function returns true if the two given float numbers are equal
bool numbers_are_equal(float a, float b, float margin)
{
	float diff=fabs(a-b);
	if (diff<=margin)
		return true;
	else
		return false;
}

// this function returns true if the S^t*S is unity and determinant(S)=+1 or -1.
bool is_S_Valid(Eigen::MatrixXd S, double delta=0.001, int LattDim=3){
	Eigen::MatrixXd SMultST(LattDim,LattDim);
	double SDet;
	SMultST=S.transpose()*S;
	SDet=S.determinant();
	//checking whether the conditions are met
	if (Matrixes_are_equal(SMultST,Unity_matrix(3),delta) && numbers_are_equal(1,fabs(SDet),delta))
		return true;
	else
		return false;
}

//this function get the matrix S and determines its proper label
std::string Label_S(Eigen::MatrixXd S,int LattDim=3)
{
	std::string label;
	double S_det=S.determinant();
	double S_trace=S.trace();
	// first we calculate the eigenvalues and then determine the number of eigenvalues equal to 1.
	Eigen::EigenSolver<Eigen::MatrixXd> ES(S,false);
	Eigen::VectorXd S_eivals = ES.eigenvalues().real();
	int count=0;
	for (int i=0;i<LattDim;i++)
		if(numbers_are_equal(S_eivals(i,0), 1, 0.001))
			count=count+1;
	//checking the conditions from which the configuration fo the given matrix will be determined
	bool det_1= numbers_are_equal(S_det, 1, 0.001);
	bool trace_3=numbers_are_equal(S_trace, 3, 0.001);
	bool trace_m3=numbers_are_equal(S_trace, -3, 0.001);
	bool det_m1= numbers_are_equal(S_det, -1, 0.001);

	if (det_1 && trace_3)
		label = "Identity\n";
	else if (det_m1 && trace_m3)
		label = "Inverse\n";
	else if (det_1 && count==1)
		label = "Rotation\n";
	else if (det_m1 && count==0 && !trace_m3)
		label = "Improper Rotation\n";
	else if (det_m1 && count==2)
		label = "Mirror\n";
	else
		printf("unrecognized configuration.\n");
	return label;

}



//this function returns true if the two given matrixes contain the same columns with any possible permutation
bool Matrixes_have_same_components(Eigen::MatrixXd latticeBasis1, Eigen::MatrixXd latticeBasis2, int NumRows, int NumCols){
	Eigen::MatrixXd sub(NumRows,1);
	int mask[NumCols];
	for (int i=0;i<NumCols;i++)
		mask[i]=1;

	for (int i=0;i<NumCols;i++)
		for(int j=0;j<NumCols;j++)
		{
			if (mask[j]==1){
			   sub=latticeBasis1.col(i)-latticeBasis2.col(j);
			   if (sub.sum()==0)
				   mask[j]=0;
			}
		}
	int mask_sum=0;
	for (int i=0;i<NumCols;i++)
		mask_sum=mask_sum+mask[i];

	if (mask_sum==0)
		return true;
	else
		return false;
}

//this function returns true f the two given vectors are identical
bool vectors_are_equal(Eigen::MatrixXd latticeBasis1, Eigen::MatrixXd latticeBasis2, int NumRows)
{
	int temp=0;
	for (int i=0;i<NumRows;i++)
		if(latticeBasis1(i)==latticeBasis2(i))
			temp=temp+1;
	//if all the elements are the same, vectors are identical.
	if (temp==NumRows)
		return true;
	else
		return false;

}

// this function returns the vector where each element is a new mesh point constructed from the given lattice
vector<vector<double>> AllPoints(Eigen::MatrixXd lattice, int radius)
{
	 vector<double> tempV{0, 0, 0};
	 vector<vector<double>> points;
	 //each new vector is the linear combination of the lattice basis vectors
	 for (int m=-radius;m<radius+1;m=m+1)
		 for (int n=-radius;n<radius+1;n=n+1)
			 for (int k=-radius;k<radius+1;k=k+1)
	 	 	 	 {
				 tempV[0]=double(m*lattice(0,0)+n*lattice(0,1)+k*lattice(0,2));
				 tempV[1]=double(m*lattice(1,0)+n*lattice(1,1)+k*lattice(1,2));
				 tempV[2]=double(m*lattice(2,0)+n*lattice(2,1)+k*lattice(2,2));
				 points.push_back(tempV);
	 	 	 	 }
     return points;
}

//this function reads the lattice or basis matrix from the file
Eigen::MatrixXd ReadingMatrixFromFile(const char * filename, int LattRow, int LattCol)
{

	Eigen::MatrixXd lattice(LattRow,LattCol);
	ifstream fp(filename);
	if (! fp) {
		cout << "Error, file couldn't be opened" << endl;
	}
	for(int row = 0; row < LattRow; row++) {  // stop loops if nothing to read
		for(int column = 0; column < LattCol; column++){
			fp >> lattice(row,column);
			if ( ! fp ) {
				cout << "Error reading file for element " << row << "," << column << endl;
			}
		}
	}
	return lattice;
}

void multiplication_table(vector<Eigen::MatrixXd> Mat)
{
	printf("\n");
	Eigen::MatrixXd multip_result;
	int mult_table[Mat.size()][Mat.size()];
	for (int i=0; i<Mat.size();i++)
		for (int j=0; j<Mat.size();j++)
		{
			multip_result=Mat[i]*Mat[j];
			for(int k=0;k<Mat.size();k++)
				if(Matrixes_are_equal(multip_result,Mat[k],0.001,3,3))
					mult_table[i][j]=k;
			//	else
			//		printf("it is not a closed group!!");
		}
	printf("\t");
	for (int i=0; i<Mat.size();i++)
			printf("%d\t",i);
	printf("\n\n");
	for (int i=0; i<Mat.size();i++)
		{
		printf("%d\t",i);
		for (int j=0; j<Mat.size();j++)
			printf("%d\t",mult_table[i][j]);
		printf("\n");
		}

}


void printing_PointGroups(vector<Eigen::MatrixXd> Mat, int LattDim=3)
{
printf("\n");
printf("************************************************************************************\n");
printf("***********************************Point Groups*************************************\n");
printf("************************************************************************************\n");
for(int i=0;i<Mat.size();i++)
	{
	cout<<Mat[i]<<endl<<endl;
	cout<<Label_S(Mat[i],LattDim)<<endl<<endl;
	printf("===========================\n");
	}
printf("************************************************************************************\n");
printf("************************************************************************************\n");
printf("************************************************************************************\n");

}

void printing_FactorGroups(vector<Eigen::MatrixXd> Mat,vector<Eigen::MatrixXd> Ta, int LattDim=3)
{
printf("\n");
printf("************************************************************************************\n");
printf("********************************Factor Groups***************************************\n");
printf("************************************************************************************\n");
for(int i=0;i<Mat.size();i++)
	{
	cout<<Mat[i]<<endl<<endl;
	cout<<Ta[i]<<endl<<endl;
	cout<<Label_S(Mat[i],LattDim)<<endl<<endl;
	printf("===========================\n");
	}
printf("************************************************************************************\n");
printf("************************************************************************************\n");
printf("************************************************************************************\n");

}



int main()
{
 //defining the used variables
	 std::string label;
	 int LattDim=3;
	 char FilePath2[200]="./Diamond_basis.txt";
	 //first checking the number of lines in the Basis file showing the number of basis
	 int numOfBasis = 0;
	 std::string line;
	 std::ifstream myfile(FilePath2);
	 while (std::getline(myfile, line))
		 ++numOfBasis;

	 Eigen::MatrixXd lattice(LattDim,LattDim);
	 Eigen::MatrixXd lattice2(LattDim,LattDim);

	 Eigen::MatrixXd latticePrime(LattDim,LattDim);
	 Eigen::MatrixXd S(LattDim,LattDim);
	 Eigen::MatrixXd latticeBasis(LattDim,numOfBasis);
	 Eigen::MatrixXd latticeBasis2(numOfBasis,LattDim);

	 Eigen::MatrixXd latticeBasisShifted(LattDim,numOfBasis);
	 Eigen::MatrixXd latticeBasisMappedS(LattDim,numOfBasis);
	 Eigen::MatrixXd latticeBasistranslation(LattDim,numOfBasis*numOfBasis);
	 vector<Eigen::MatrixXd> PointGroup;
	 vector<Eigen::MatrixXd> FactorGroup_s;
	 vector<Eigen::MatrixXd> FactorGroup_ta;


	 //reading the lattice from the file
	 char FilePath[200]="./Diamond_lattice.txt";
	 lattice2=ReadingMatrixFromFile(FilePath, LattDim, LattDim);
	 latticeBasis2=ReadingMatrixFromFile(FilePath2, numOfBasis,LattDim);
	 lattice=lattice2.transpose();
	 latticeBasis=latticeBasis2.transpose();

	 cout<<"given lattice is:"<<endl<<lattice<<endl<<endl;
	 cout<<"given basis is:"<<endl<<latticeBasis<<endl<<endl;
	 cout<<"==========================================================\n";

	 double delta=0.001;
	 int radius=2;
	 //make the mesh from the lattice
	 vector<vector<double>> points=AllPoints(lattice, radius);
	 int numberOfValidS=0;
	 int numberOfValidSandT=0;
	 int temp=0;
	 int index1=0;

	 //creating lattice primes, then constructing S, and checking if it is a valid S.
	 for (int i=0; i<points.size();i=i+1)
		 for (int j=0; j<points.size();j=j+1)
			 for (int k=0; k<points.size();k=k+1)

		    {
	 	 	 latticePrime(0,0)=points[i][0];
	 	 	 latticePrime(1,0)=points[i][1];
	 	 	 latticePrime(2,0)=points[i][2];

	 	 	 latticePrime(0,1)=points[j][0];
	 	 	 latticePrime(1,1)=points[j][1];
	 	 	 latticePrime(2,1)=points[j][2];

	 	 	 latticePrime(0,2)=points[k][0];
	 	 	 latticePrime(1,2)=points[k][1];
	 	 	 latticePrime(2,2)=points[k][2];


	 	 	 S=latticePrime*lattice.inverse();
	 	 	 if(is_S_Valid(S,delta,LattDim))
	 	 	 	 {
	 	 		 PointGroup.push_back(S);
	 	 		 numberOfValidS=numberOfValidS+1;
	 	 		 //now that S is valid, we can check the basis
	 	 		 latticeBasisMappedS=S*latticeBasis;
	 	 		 index1=0;
	 	 		 //calculating all the possible translations
	 	 		 for (int i2=0; i2<numOfBasis;i2=i2+1)
	 	 			 for (int j2=0; j2<numOfBasis;j2=j2+1)
	 	 			 {
	 	 				 latticeBasistranslation.col(index1)=latticeBasisMappedS.col(i2)-latticeBasis.col(j2);
	 	 				 index1=index1+1;
	 	 			 }
	 	 		//for each of the possible translation, we check to see whether r and r' can be mapped together
	 	 		for (int i2=0; i2<numOfBasis*numOfBasis; i2=i2+1)
	 	 		{
	 	 			//because we can have several similar translations, we exclude the duplicates fisrt
	 	 			temp=0;
	 	 			for (int i3=0;i3<i2;i3++)
	 	 				if(vectors_are_equal(latticeBasistranslation.col(i2),latticeBasistranslation.col(i3),LattDim))
	 	 					temp=1;
	 	 			if (temp==0){
	 	 				//if the translation is not a duplicate, we first shift the original basis with the translation and check whether it mapps to the maped basis.
	 	 				latticeBasisShifted=latticeBasis.colwise()+latticeBasistranslation.col(i2);
	 	 				//we also exclude the zero translation.
	 	 				//if((latticeBasistranslation.col(i2).sum()!=0)&&Matrixes_have_same_components(latticeBasisMappedS,latticeBasisShifted,LattDim,numOfBasis))
	 	 				if(Matrixes_have_same_components(latticeBasisMappedS,latticeBasisShifted,LattDim,numOfBasis))
	 	 					{numberOfValidSandT=numberOfValidSandT+1;
	 	 					//pushing back the facto groups
	 	 					FactorGroup_s.push_back(S);
	 	 					FactorGroup_ta.push_back(latticeBasistranslation.col(i2));
	 	 					//printing Ta
	 	 					//cout<<latticeBasistranslation.col(i2)<<endl<<endl;
	 	 					//printing S
	 	 					//cout<<S<<endl<<endl;
	 	 					//labeling the S
	 	 					//label=Label_S(S,LattDim);
	 	 					//cout<<label<<endl;
	 	 					//cout<<"=========================================================="<<endl;
	 	 					}
	 	 						}
	 	 		}
	 	 	 	}//is_S_Valid

		    }
	 printing_PointGroups(PointGroup,LattDim);
	 printing_FactorGroups(FactorGroup_s,FactorGroup_ta,LattDim);
	 cout<<"Total Number of Valid Point Groups:  "<<numberOfValidS<<endl;
	 cout<<"Total number of valid Factor Groups:   "<<numberOfValidSandT<<endl;
	 multiplication_table(PointGroup);


return 0;
}
