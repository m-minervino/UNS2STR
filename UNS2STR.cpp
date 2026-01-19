#include <iostream>
#include <fstream>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

//Required by TecIO
#include "/opt/tecplot/360ex_2025r2/include/TECIO.h"
#include <cstring>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <stdexcept>
#include <vector>

using namespace std;
int main(int argc, char** argv) {
	void* inputFileHandle = NULL;
	//Set debug mode FLAG
	    bool debug_mode = 0;
	    if (argc <= 3) {
	    	debug_mode = false;
	    }
	    else if (atof(argv[3]) == 1) {
	        debug_mode = true;
	    }
    //Set z-overx scaling factor
	    double scale_factor = 1;
	    if (argc > 4) {	    	
	        scale_factor = atof(argv[4]);
            cout << "scale factor:" << scale_factor << "\n";
	    }
    //Set grid motion
        /* double x_orig=0;    double z_orig=0;
        double ampl_x=0;    double freq_x=0;    double t0_x=0;
        double ampl_z=0;    double freq_z=0;    double t0_z=0;
        double ampl_a=0;    double freq_a=0;    double t0_a=0;
        double dt=0;        int t_step=0;
        ifstream motion_file(argv[3]);
        if (motion_file.is_open()) {
            string line;
            int i = 0;
            int m = 0;
            while (getline(motion_file, line)) {
                i = line.find_last_of("\t ");
                m = line.find_first_of("\t ");
                if (line.substr(0, m) == "X_ORIG") {
                    x_orig = atof(line.substr(i + 1).data());
                }
                else if (line.substr(0, m) == "Z_ORIG") {
                    z_orig = atof(line.substr(i + 1).data());
                }
                else if (line.substr(0, m) == "AMPL_X") {
                    ampl_x = atof(line.substr(i + 1).data());
                }
                else if (line.substr(0, m) == "FREQ_X") {
                    freq_x = atof(line.substr(i + 1).data());
                }
                else if (line.substr(0, m) == "T0_X") {
                    t0_x = atof(line.substr(i + 1).data());
                }
                else if (line.substr(0, m) == "AMPL_Z") {
                    ampl_z = atof(line.substr(i + 1).data());
                }
                else if (line.substr(0, m) == "FREQ_Z") {
                    freq_z = atof(line.substr(i + 1).data());
                }
                else if (line.substr(0, m) == "T0_Z") {
                    t0_z = atof(line.substr(i + 1).data());
                }
                else if (line.substr(0, m) == "AMPL_A") {
                    ampl_a = atof(line.substr(i + 1).data());
                }
                else if (line.substr(0, m) == "FREQ_A") {
                    freq_a = atof(line.substr(i + 1).data());
                }
                else if (line.substr(0, m) == "T0_A") {
                    t0_a = atof(line.substr(i + 1).data());
                }
                else if (line.substr(0, m) == "DT") {
                    dt = atof(line.substr(i + 1).data());
                }
                else if (line.substr(0, m) == "T_STEP") {
                    t_step = atof(line.substr(i + 1).data());
                }
            }
            motion_file.close();
        } */
	//Read structured grid
	int R = tecFileReaderOpen((char const*)(argv[1]), &inputFileHandle);
    int64_t iMax=0;
	int64_t kMax=0;
    int64_t tmp_max = 0;
    R = tecZoneGetIJK(inputFileHandle, 1, &iMax, &tmp_max, &kMax);
	if (kMax==1) {
		R = tecZoneGetIJK(inputFileHandle, 1, &iMax, &kMax, &tmp_max);
	}
    if (debug_mode) {
        cout << "iMax is:   " << iMax << "\n";
        cout << "kMax is:   " << kMax << "\n";
    }
    int32_t numVars = 0;
    R = tecDataSetGetNumVars(inputFileHandle, &numVars);
    vector < int64_t > numValues(numVars,0);
    vector < string > in_var_names(numVars,"");
    for (int j = 0; j < numVars; j++) {
        char* name = NULL;
        R = tecVarGetName(inputFileHandle, j + 1, &name);
        in_var_names[j] = name;
    }
    int64_t max_nv = 0;
    for (int j = 0; j < numVars; j++) {
        R = tecZoneVarGetNumValues(inputFileHandle, 1, j + 1, &numValues[j]);
        if (numValues[j] > max_nv) {
            max_nv = numValues[j];
        }
    }
    if (debug_mode) {
        for (int j = 0; j < numVars; j++) {
            cout << "Values of zone 1 for variable " << j << " is:   " << numValues[j] << "\n";
        }
        cout << "Maximum number of variable values is:   " << max_nv << "\n";
    }
    vector < vector < double > > values(numVars,vector < double > (max_nv,0));
    for (int j = 0; j < numVars; j++) {
        R = tecZoneVarGetDoubleValues(inputFileHandle, 1, j + 1, 1, numValues[j], &values[j][0]);
    }
    R = tecFileReaderClose(&inputFileHandle);

    //Reshape data based on a i-k-var 3D matrix
    vector < vector < vector < double > > > STR_values(iMax, vector < vector < double > > (kMax, vector < double > (4,0)));
    for (int j = 0; j < numVars; j++) {
        for (int k = 0; k < kMax; k++) {
            for (int m = 0; m < iMax; m++) {
                STR_values[m][k][j] = values[j][k * iMax + m];
            }
        }
    }
    /* if (numVars==3) {   //No topologycal information is provided. Need to match nodes...
        //Apply grid motion at selected time step
        double AoA=0;
        AoA=(M_PI/180)*ampl_a*sin(freq_a*(t_step*dt+t0_a));
        for (int k = 0; k < kMax; k++) {
            for (int m = 0; m < iMax; m++) {
                STR_values[m][k][0] = ((STR_values[m][k][0]-x_orig)*cos(AoA)-(STR_values[m][k][2]-z_orig)*sin(AoA))+ampl_x*sin(freq_x*(t_step*dt+t0_x));
                STR_values[m][k][2] = ((STR_values[m][k][0]-x_orig)*sin(AoA)+(STR_values[m][k][2]-z_orig)*cos(AoA))+ampl_z*sin(freq_z*(t_step*dt+t0_z));
            }
        }
    } */

	//Read SU2 solution
	R = tecFileReaderOpen((char const*)(argv[2]), &inputFileHandle);
    int32_t numVars_SU2 = 0;
    R = tecDataSetGetNumVars(inputFileHandle, &numVars_SU2);
    vector < int64_t > numValues_SU2(numVars_SU2,0);
    vector < string > in_var_names_SU2(numVars_SU2,"");
    for (int j = 0; j < numVars_SU2; j++) {
        char* name = NULL;
        R = tecVarGetName(inputFileHandle, j + 1, &name);
        in_var_names_SU2[j] = name;
    }
    int64_t max_nv_SU2 = 0;
    for (int j = 0; j < numVars_SU2; j++) {
        R = tecZoneVarGetNumValues(inputFileHandle, 1, j + 1, &numValues_SU2[j]);
        if (numValues_SU2[j] > max_nv_SU2) {
            max_nv_SU2 = numValues_SU2[j];
        }
    }
    if (debug_mode) {
        for (int j = 0; j < numVars_SU2; j++) {
            cout << "Values of zone 1 for SU2 variable " << j << " is:   " << numValues_SU2[j] << "\n";
        }
        cout << "Maximum number of SU2 variable values is:   " << max_nv_SU2 << "\n";
    }
    vector < vector < double > > values_SU2(numVars_SU2,vector < double > (max_nv_SU2,NAN));
    for (int j = 0; j < numVars_SU2; j++) {
        R = tecZoneVarGetDoubleValues(inputFileHandle, 1, j + 1, 1, numValues_SU2[j], &values_SU2[j][0]);
    }
    double time=NAN;
    R = tecZoneGetSolutionTime(inputFileHandle, 1, &time);
    R = tecFileReaderClose(&inputFileHandle);

	//Fill-in new array of ordered unstructured data
    vector < vector < vector < double > > > ORD_values(iMax, vector < vector < double > > (kMax, vector < double > (numVars_SU2,NAN)));
	int rho_id=-1;
	for (int j = 0; j < numVars_SU2; j++) {
		if (in_var_names_SU2[j]=="Density") {
			rho_id=j;
		}
	}
	if (debug_mode) {
		cout << "Variable ID for Density is: " << rho_id << "\n";
	}
    for (int k = 0; k < kMax; k++) {
        for (int m = 0; m < iMax; m++) {
            if (numVars==3) {   //No topologycal information is provided. Need to match nodes...
                double min_dist=100000000;
                double dist=100000000;
                STR_values[m][k][3]=NAN;
                for (int i=0; i<max_nv_SU2; i++) {                
                    dist=sqrt(pow(values_SU2[0][i]-STR_values[m][k][0],2)+pow(scale_factor*(values_SU2[2][i]-STR_values[m][k][2]),2));
                    if (dist<min_dist) {
                        min_dist=dist;
                        STR_values[m][k][3]=i;
                    }
                }
                if (debug_mode) {
                    cout << "Identified node: " << k*iMax+m+1 << " of " << max_nv << "\n";
                }
            }
            for (int j = 0; j < numVars_SU2; j++) {
                if ((in_var_names_SU2[j]=="Momentum_x")||(in_var_names_SU2[j]=="Momentum_y")||(in_var_names_SU2[j]=="Momentum_z")) {
                    ORD_values[m][k][j] = values_SU2[j][int(STR_values[m][k][3])]/values_SU2[rho_id][int(STR_values[m][k][3])];
                } else {
                    ORD_values[m][k][j] = values_SU2[j][int(STR_values[m][k][3])];
                }
            }
        }
    }
	for (int j = 0; j < numVars_SU2; j++) {
        if (in_var_names_SU2[j]=="Momentum_x") {
            in_var_names_SU2[j]="Velocity_x";
		} else if (in_var_names_SU2[j]=="Momentum_y") {
            in_var_names_SU2[j]="Velocity_y";
        } else if (in_var_names_SU2[j]=="Momentum_z") {
            in_var_names_SU2[j]="Velocity_z";
        }
    }
	if (debug_mode) {
		for (int k = 0; k < kMax; k++) {
	        for (int m = 0; m < iMax; m++) {
				if (isnan(ORD_values[m][k][0])) {
					cout << "NAN detected at position:   " << m << " , " << k << "\n";
				}
			}
		}
	}
    //SU2 output data structure
    vector < vector < double > > values_SU2_NEW(numVars_SU2,vector < double > (iMax*kMax,0));
	for (int j = 0; j < numVars_SU2; j++) {
		for (int k = 0; k < kMax; k++) {
    	    for (int m = 0; m < iMax; m++) {
				values_SU2_NEW[j][iMax*k+m]=ORD_values[m][k][j];
			}
		}
	}
	//Write SU2 ordered Tecplot file
	void* outputFileHandle = NULL;
	int32_t outputDebugInfo = 1;
	ostringstream outputStream;
        for (int32_t var = 0; var < numVars_SU2; var++)
        {
            outputStream << in_var_names_SU2[var];
            if (var < numVars_SU2-1) {
                outputStream << ',';
			}
        }
	int32_t fileFormat = 1; // .szplt
    vector < int32_t > varTypes(numVars_SU2,2);
    vector < int32_t > shareVarFromZone(numVars_SU2,0);
	vector < int32_t > valueLocation(numVars_SU2,1);
	vector < int32_t > passiveVarList(numVars_SU2,0);
	int32_t shareConnectivityFromZone=0;
	int64_t numFaceConnections=0;
	int32_t faceNeighborMode=0;
	int32_t outputZone;
	R = tecFileWriterOpen("SU2.szplt", "Ordered_SU2_solution", outputStream.str().c_str(), fileFormat, 0, 1, NULL, &outputFileHandle);
	R = tecFileSetDiagnosticsLevel(outputFileHandle, outputDebugInfo);
	R = tecZoneCreateIJK(outputFileHandle, "Ordered_SU2_data", iMax, 1, kMax, &varTypes[0],
                    &shareVarFromZone[0], &valueLocation[0], &passiveVarList[0],
                    shareConnectivityFromZone, numFaceConnections, faceNeighborMode, &outputZone);
	for (int j = 0; j < numVars_SU2; j++) { 
		R = tecZoneVarWriteDoubleValues(outputFileHandle, outputZone, j+1, 0, iMax*kMax, &values_SU2_NEW[j][0]);
	}
    R = tecZoneSetUnsteadyOptions(outputFileHandle, outputZone, time, 1);
	R = tecFileWriterClose(&outputFileHandle);

    //Conditionally write-out input structured grid with topologycal information included
    if (numVars==3) {
        //Output data structure
        vector < vector < double > > values_NEW(4,vector < double > (iMax*kMax,0));
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < kMax; k++) {
                for (int m = 0; m < iMax; m++) {
                    values_NEW[j][iMax*k+m]=values[j][k * iMax + m];
                }
            }
        }
        for (int k = 0; k < kMax; k++) {
            for (int m = 0; m < iMax; m++) {
                values_NEW[3][iMax*k+m]=STR_values[m][k][3];
            }
        }
        //Write output data
        ostringstream outputStream_STR;
            for (int32_t var = 0; var < 3; var++) {
                outputStream_STR << in_var_names[var];
                outputStream_STR << ',';
            }
            outputStream_STR << "SU2_NODE-ID";
        int32_t varTypes_STR[4];
        for (int i=0; i<4; i++) {
            varTypes_STR[i]=2;
        }
        int32_t shareVarFromZone_STR[4];
        for (int i=0; i<4; i++) {
            shareVarFromZone_STR[i]=0;
        }
        int32_t valueLocation_STR[4];
        for (int i=0; i<4; i++) {
            valueLocation_STR[i]=1;
        }
        int32_t passiveVarList_STR[4];
        for (int i=0; i<4; i++) {
            passiveVarList_STR[i]=0;
        }
        R = tecFileWriterOpen("TOPO.szplt", "Structured grid with SU2 topo map", outputStream_STR.str().c_str(), fileFormat, 0, 1, NULL, &outputFileHandle);
        R = tecFileSetDiagnosticsLevel(outputFileHandle, outputDebugInfo);
        R = tecZoneCreateIJK(outputFileHandle, "DATA", iMax, 1, kMax, &varTypes_STR[0],
                        &shareVarFromZone_STR[0], &valueLocation_STR[0], &passiveVarList_STR[0],
                        shareConnectivityFromZone, numFaceConnections, faceNeighborMode, &outputZone);
        for (int j = 0; j < 4; j++) { 
            R = tecZoneVarWriteDoubleValues(outputFileHandle, outputZone, j+1, 0, iMax*kMax, &values_NEW[j][0]);
        }
        R = tecFileWriterClose(&outputFileHandle);
    }
}
