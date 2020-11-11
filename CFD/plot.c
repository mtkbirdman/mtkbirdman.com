#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define GNUPLOT_PATH "C:/gnuplot/bin/gnuplot.exe" //PATH for gnuplot.exe

int main()
{
	FILE *gp,*fp; //gp:gnuplot pointer, fp:file pointer
	char foil_name[256],file_name[256],foil_data[256],file_format[]=".txt",*ch,str[256];
	int i,j,file_size=0;
	double x,y,val[1024][4];

	//Input foil name
	printf("Input foil name\n>>");
	fgets(foil_name,sizeof(foil_name),stdin); //Input from keybord
	foil_name[strlen(foil_name)-1]='\0'; //Remove "\n"
	sprintf(file_name,"%s%s",foil_name,file_format); //Concatinate foil_name & file_format

	//Open file
	if((fp=fopen(file_name,"r"))==NULL){ //Failed
		fprintf(stderr,"Not Found %s.",file_name); //Standard error output
		exit(EXIT_FAILURE); //Program abnormal termination
	}else{ //Success
		while(fgets(foil_data,sizeof(foil_data),fp)!=NULL){ //Read the data to the last row
			ch=strtok(foil_data," \n"); //Separate the first value by a delimiter
			for(i=0;i<256;i++){
				if(ch==NULL){ //End of string
					break;
				}else{
					val[file_size][i]=atof(ch); //Convert a read value from char to double
				}
				ch=strtok(NULL," \n"); //Separate the next value
			}
			file_size++; //Update file size
		}
	}
	fclose(fp); //Close file

	//Start Gnuplot comand
	if((gp=_popen(GNUPLOT_PATH,"w"))==NULL){ //Start gnuplot with pipe
		fprintf(stderr,"Not Found %s.",GNUPLOT_PATH); //Standard error output
		exit(EXIT_FAILURE); //Program abnormal termination
	}

	//Send comands to gnuplot
	fprintf(gp, "set xrange [0:1]\n"); //Set ranges
	fprintf(gp, "set yrange [-0.5:0.5]\n");
	fprintf(gp,"set size square\n");
	fprintf(gp, "plot '-' with lines\n"); //Plot array

	for(i=1;i<=file_size-1;i++){
	   fprintf(gp,"%f\t%f\n",val[i][0],val[i][1]); //Write data
    }
    fprintf(gp,"e\n"); //End of array

	fflush(gp); //Spit out the data stored in the buffer (required)
	system("pause");
	fprintf(gp, "exit\n"); //Terminate gnuplot
	_pclose(gp); //Close Gnuplot
	
}