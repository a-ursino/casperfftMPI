#include <iostream>
#include <cmath>
#include "cooleyTukey.h"
#include "fftutils.h"

using namespace std;

/* Functions Implementation*/

void transpose2(double *data3DRr,double *data3DRi, int xR, int yR, int off) 
{
	int xxx, var1 = 0, var2 = 0;
	double temp1,temp2;

	for (xxx = 0; xxx <= xR * yR - 1; xxx++) {	// for each position
		var1 = xxx;
		var2 = 0;
		do {	
			var2++;
			var1 = (var1 % yR) * xR + var1 / yR;
		} while (var1 > xxx);
		if (var1 < xxx || var2 == 1) continue;
		var1 = xxx;	
		temp1 = data3DRr[(off+var1)];
		temp2 = data3DRi[(off+var1)];
		do {
			var2 = (var1 % yR) * xR + var1 / yR;
			data3DRr[(off+var1)]   = (var2 == xxx) ? temp1 : data3DRr[(off+var2)];
			data3DRi[(off+var1)] = (var2 == xxx) ? temp2 : data3DRi[(off+var2)];
			var1 = var2;
		} while (var1 > xxx);
	}
}

void transpose3(const unsigned offset, double *data3DRr, double *data3DRi, int zR, int yR, int xR, int off) 
{
	double *slice = (double *) malloc(sizeof(double)*yR*zR);
	double *slice1 = (double *) malloc(sizeof(double)*yR*zR);
	int counter;
	int newX, newY, newZ;
	int nX = xR, nY = zR;//, nZ = yR;
	for (newZ = 0; newZ < xR; newZ++) { // yR
		counter = 0;
		for (newY = 0; newY < yR; newY++) {
			for (newX = zR - 1; newX >= 0; newX--) {
				slice[counter] = data3DRr[(offset + off + newZ + newY * xR + newX * xR * yR)];
				slice1[counter] = data3DRi[(offset + off + newZ + newY * xR + newX * xR * yR)];
				counter++;
			}
		}
		transpose2(slice,slice1, zR, yR, 0);
		counter = 0;
		for (newY = zR - 1; newY >= 0; newY--) {
			for (newX = 0; newX < yR; newX++) {
				//data3DRr[(offset + off + newZ + newY * xR + newX * xR * yR)] = slice[counter];
				//data3DRi[(offset + off + newZ + newY * xR + newX * xR * yR)] = slice1[counter];
				data3DRr[(offset + off + newZ + newY * nX + newX * nX * nY)] = slice[counter];
				data3DRi[(offset + off + newZ + newY * nX + newX * nX * nY)] = slice1[counter];
				counter++;
			}
		}
	}
	free(slice);
	free(slice1);
}
void convolveCPU(const unsigned offset1,int ASPAN){
	//if (useCpu == 0) return; 
	for (int iz = 0; iz < xRange * yRange * zRange; iz++) {

		hrhVecI[(iz+offset1)] = hrRaVec[(0*ASPAN+iz+offset1)]   * hrRmVecI[iz+offset1] - hiRaVec[(0*ASPAN+iz+offset1)] 	* hiRmVecI[iz+offset1] +	// XX
								hrRaVec[(1*ASPAN+iz+offset1)] 	* hrRmVecJ[iz+offset1] - hiRaVec[(1*ASPAN+iz+offset1)] 	* hiRmVecJ[iz+offset1] +	// XY
								hrRaVec[(2*ASPAN+iz+offset1)]   * hrRmVecK[iz+offset1] - hiRaVec[(2*ASPAN+iz+offset1)] 	* hiRmVecK[iz+offset1] ;	// XZ

		hihVecI[(iz+offset1)] = hiRaVec[(0*ASPAN+iz+offset1)] 	* hrRmVecI[iz+offset1] + hrRaVec[(0*ASPAN+iz+offset1)]   * hiRmVecI[iz+offset1] +	// XX
								hiRaVec[(1*ASPAN+iz+offset1)] 	* hrRmVecJ[iz+offset1] + hrRaVec[(1*ASPAN+iz+offset1)]   * hiRmVecJ[iz+offset1] +	// XY
								hiRaVec[(2*ASPAN+iz+offset1)] 	* hrRmVecK[iz+offset1] + hrRaVec[(2*ASPAN+iz+offset1)]   * hiRmVecK[iz+offset1];	// XZ

		hrhVecJ[(iz+offset1)] = hrRaVec[(1*ASPAN+iz+offset1)]   * hrRmVecI[iz+offset1] - hiRaVec[(1*ASPAN+iz+offset1)] 	* hiRmVecI[iz+offset1] +	// YX
								hrRaVec[(3*ASPAN+iz+offset1)] 	* hrRmVecJ[iz+offset1] - hiRaVec[(3*ASPAN+iz+offset1)] 	* hiRmVecJ[iz+offset1] +	// YY
								hrRaVec[(4*ASPAN+iz+offset1)]   * hrRmVecK[iz+offset1] - hiRaVec[(4*ASPAN+iz+offset1)] 	* hiRmVecK[iz+offset1] ;	// YZ

		hihVecJ[(iz+offset1)] = hiRaVec[(1*ASPAN+iz+offset1)] 	* hrRmVecI[iz+offset1] + hrRaVec[(1*ASPAN+iz+offset1)]   * hiRmVecI[iz+offset1] +	// YX
								hiRaVec[(3*ASPAN+iz+offset1)] 	* hrRmVecJ[iz+offset1] + hrRaVec[(3*ASPAN+iz+offset1)]   * hiRmVecJ[iz+offset1] +	// YY
								hiRaVec[(4*ASPAN+iz+offset1)] 	* hrRmVecK[iz+offset1] + hrRaVec[(4*ASPAN+iz+offset1)]   * hiRmVecK[iz+offset1];	// YZ

		hrhVecK[(iz+offset1)] = hrRaVec[(2*ASPAN+iz+offset1)]   * hrRmVecI[iz] - hiRaVec[(2*ASPAN+iz+offset1)] 	* hiRmVecI[iz+offset1] +			// ZX
								hrRaVec[(4*ASPAN+iz+offset1)] 	* hrRmVecJ[iz+offset1] - hiRaVec[(4*ASPAN+iz+offset1)] 	* hiRmVecJ[iz+offset1] +	// ZY
								hrRaVec[(5*ASPAN+iz+offset1)]   * hrRmVecK[iz+offset1] - hiRaVec[(5*ASPAN+iz+offset1)] 	* hiRmVecK[iz+offset1] ;	// ZZ

		hihVecK[(iz+offset1)] = hiRaVec[(2*ASPAN+iz+offset1)] 	* hrRmVecI[iz+offset1] + hrRaVec[(2*ASPAN+iz+offset1)]   * hiRmVecI[iz+offset1] +	// ZX
								hiRaVec[(4*ASPAN+iz+offset1)] 	* hrRmVecJ[iz+offset1] + hrRaVec[(4*ASPAN+iz+offset1)]   * hiRmVecJ[iz+offset1] +	// ZY
								hiRaVec[(5*ASPAN+iz+offset1)] 	* hrRmVecK[iz+offset1] + hrRaVec[(5*ASPAN+iz+offset1)]   * hiRmVecK[iz+offset1];	// ZZ
	}
	int show_results=1;
	if (show_results){
		printf("CONVOLUTION hVecI\n");
		printf("XRange=:%d\t YRange=:%d\t ZRange=:%d\t  \n", xRange,yRange,zRange);
		printMe(0,hrhVecI, hihVecI, zRange, yRange, xRange, 0 );
		/*printf("CONVOLUTION hVecJ\n");
		printMe(0,hrhVecJ, hihVecJ, zRange, yRange, xRange, 0 );
		printf("CONVOLUTION hVecK\n");
		printMe(0,hrhVecK, hihVecK, zRange, yRange, xRange, 0 );*/
		show_results=0;
	}
}

void cooleyTukeyCpu3DFFT(const unsigned offset, const unsigned  N, const unsigned size,double *data3DFr,double *data3DFi,double *data3DRr,double *data3DRi,int ASPAN_offset, int show_results, int fft_type){
	if (size == 0) return;
	//if (useCpu == 0) return;
	int planeStart=0;
	int tHolder = 0;

	//reinitializing them when the ranges are different
	if(fft_type !=1){
		xRange = initialxRange;
		yRange = initialyRange;
		zRange = initialzRange;
	}
	printf("XRange=:%d\t YRange=:%d\t ZRange=:%d\t  \n", xRange,yRange,zRange);
	if (show_results == 1) {
		printf("Start\n");
		printMe(offset,data3DFr, data3DFi, zRange, yRange, xRange, ASPAN_offset );
	}
		
	//************************** X TRANSFORM STARTS ********************************************//
	cooleyTukeyCpu(offset, xRange, size,data3DFr,data3DFi,data3DRr,data3DRi,ASPAN_offset,show_results,fft_type,0,0,0);

	if (show_results)
	{
		printf("X-TRANSFORM \n");
		printMe(offset,data3DRr, data3DRi, zRange, yRange, xRange, ASPAN_offset );
		//show_results=0;
	}
	//************************** X TRANSFORM ENDS ********************************************//

	//************************** XY PLANE TRANSPOSE STARTS *******************************************//
	if (yRange > 1)
	{
		for (int z = 0; z < zRange; z++) 
		{
			planeStart = offset + ASPAN_offset + z * xRange * yRange;
			transpose2(data3DRr,data3DRi, xRange, yRange, planeStart);
		}


		if (show_results)
		{
			printf("Transpose XY \n");
			printMe(offset,data3DRr, data3DRi, zRange, yRange, xRange, ASPAN_offset );
		//show_results=0;
		}

	}
	//************************** XY PLANE TRANSPOSE ENDS *******************************************//

	//************************** Y TRANSFORM STARTS ********************************************//
	if (yRange > 1)
	{
		// N=yRange; //xRange=yRange;
		// and yRange=xRange;
		tHolder = xRange; xRange = yRange; yRange = tHolder;

		printf("Y-TRANSFORM \n");
		printf("XRange=:%d\t YRange=:%d\t ZRange=:%d\t  \n", xRange,yRange,zRange);

		cooleyTukeyCpu(offset, xRange, size,data3DRr,data3DRi,data3DFr,data3DFi,ASPAN_offset,show_results,fft_type,0,0,0);

		if (show_results)
		{
			printf("Y-TRANSFORM \n");
			printMe(offset,data3DFr, data3DFi, zRange, yRange, xRange, ASPAN_offset );
		//show_results=0;
		}
	}
	//************************** Y TRANSFORM ENDS ********************************************//

	//************************** YZ PLANE TRANSPOSE STARTS *******************************************//
	if (zRange > 1)	
	{
		transpose3(offset,data3DFr, data3DFi, zRange, yRange, xRange, ASPAN_offset);
		tHolder = zRange; zRange = yRange; yRange = tHolder;
		printf("Transpose YZ \n");
		printf("XRange=:%d\t YRange=:%d\t ZRange=:%d\t  \n", xRange,yRange,zRange);
		if (show_results)
		{
			printf("Transpose YZ \n");
			printMe(offset,data3DFr, data3DFi, zRange, yRange, xRange, ASPAN_offset );
		//show_results=0;
		}
	}
	//************************** YZ PLANE TRANSPOSE ENDS *******************************************//

	//************************** ZX PLANE TRANSPOSE STARTS *******************************************//
	if (zRange > 1)	
	{
		for (int z = 0; z < zRange; z++)
		{
			planeStart =  offset + ASPAN_offset + z * xRange * yRange;
			transpose2(data3DFr, data3DFi, xRange, yRange, planeStart);
		}
		tHolder = xRange; xRange = yRange; yRange = tHolder;
		printf("Transpose ZX \n");
		printf("XRange=:%d\t YRange=:%d\t ZRange=:%d\t  \n", xRange,yRange,zRange);
		if (show_results)
		{
			printf("Transpose ZX \n");
			printMe(offset,data3DFr, data3DFi, zRange, yRange, xRange, ASPAN_offset );
		//show_results=0;
		}
	}
	//************************** ZX PLANE TRANSPOSE ENDS *******************************************//

	//************************** Z TRANSFORM STARTS ********************************************//
	if (zRange > 1)	
	{
		cooleyTukeyCpu(offset, xRange, size,data3DFr,data3DFi,data3DRr,data3DRi,ASPAN_offset,show_results,fft_type,0,0,0);

		// normalization in inverse fft starts

		if (fft_type == 1) {
			for (int z = 0; z < zRange; z++) {
				int planeStart = z * xRange * yRange;
				for (int y = yRange-1; y >= 0; y--) {
					int yStart = y * xRange;
					for (int x = 0; x < xRange; x++) {
						data3DRr[(offset+ASPAN_offset+planeStart+yStart+x)] = data3DRr[(offset+ASPAN_offset+planeStart+yStart+x)]/(xRange*yRange*zRange);
						data3DRi[(offset+ASPAN_offset+planeStart+yStart+x)] = data3DRi[(offset+ASPAN_offset+planeStart+yStart+x)]/(xRange*yRange*zRange);
					}
				}
			}

		}
		// normalization in inverse fft ends

		if (show_results)
		{
			printf("Z-TRANSFORM \n");
			printMe(offset,data3DRr, data3DRi, zRange, yRange, xRange, ASPAN_offset );
			//system("rm plotMe && touch plotMe");
			//show_results=0;
		}
	}
	//************************** Z TRANSFORM ENDS ********************************************//
	//	printField(hVecI, hVecJ, hVecK, zRange, yRange, xRange);
}

//   Inplace version of rearrange function
void cooleyTukeyCpu(const unsigned offset, const unsigned  N, const unsigned size,double *data3DFr,double *data3DFi,double *data3DRr,double *data3DRi,int ASPAN_offset, int show_results, int fft_type, int full, int red, int reduction)
{
	//printf("red = %d \n", red); 
	//printf("xRange=%d \n",N);
	if (size == 0) return;
	//if (useCpu == 0) return; 
	const unsigned powN = (unsigned)log2(N);
	//const double start = omp_get_wtime();
	int red_off=0;
	//TODO:: set the number of threads
	// #pragma omp parallel for
	//*********************************3dFFT bit revrsal for N or aVec Starts*******************************
	// start of bit reversal see chapter printed page No. 465
	int red2=red*N; 
	for (int i = 0; i < (int)size; ++i) {

   		// for reduction
		if(i !=0 && reduction==1 && i % (red*N) ==0 )
		{
			if(red % 2==0 )
				i =i + ((full-red)*N);
			else if(i % red2 ==0 && red % 2!=0)
			{
				i =i + ((full-red)*N);
				red2+=(red*N)+((full-red)*N);
			}

		}
		unsigned int lIndex =  i % N;
		unsigned int lPosition  = 0;
		unsigned int lReverse= 0;
		while(lIndex) {
			lReverse = lReverse << 1;
			lReverse += lIndex %2;
			lIndex = lIndex>>1;
			lPosition++;
		}
		if (lPosition < powN) {
            lReverse = lReverse << (powN-  lPosition);   //<< is a shift operator
        }
        uint lTempReverse = lReverse + (i / N) * N + offset;

        data3DRr[lTempReverse + ASPAN_offset] = data3DFr[i + offset + ASPAN_offset];
        data3DRi[lTempReverse + ASPAN_offset] = data3DFi[i + offset + ASPAN_offset];
    }

	//*********************************3dFFT bit revrsal for N or aVec Ends*******************************

	/* Some Complex Arithmetic for reference
	 * (a+bi) + (c+di) = (a+c) + (b+d)i
	 * (a+bi) - (c+di) = (a-c) + (b-d)i
	 * (a+bi) * (c+di) = (ac - bd) + (bc + ad)i
	 * (a+bi) / (c+di) = [ (ac + bd) / (c^2 + d^2) ] + [ (bc - ad) / (c^2 + d^2) ]i
	 */

	//double pi = acos(-1);		// value of pi .. never use 22/7 again
	const double twopi =  2 * 3.14159265358979323846;
	double cs;
	double sn;
    double norm = 1.0/(xRange*yRange*zRange);  //N;   // 'normalisation' factor (xRange*yRange*zRange);
    int red1=red;             
    red_off=0;
    for (unsigned i = 0; i < size / N; ++i) {        //Number of FFTS      //        In kernel it is = unsigned int lIndex =  lThread %n;

   		//  printf("red = %d \n", red); 

   		// for reduction
    	if(i !=0 && reduction==1 &&  i % red==0 )
    	{
    		if(red % 2==0){
    			i = i + (full-red);
    		}
    		else if(red % 2!=0 && i==red){
				i = i + (full-red);
    			red1+=red+(full-red);
    		}
        	// if(i>size / N)
        	// return;
    	}
    	if(i !=0 && reduction==1 &&  i  % red1==0 && red % 2!=0 && (i-1) > red)
    	{
			i = i + (full-red);//(full)-(red);
			red1+=red+(full-red);
			// if(i>size / N)
			// return;
		}



        for (unsigned p = 0; p < powN ; ++p ) {            // for number of stages         // In kernel it is = for(Iter =0 ,nIter =1;Iter<(powN);Iter ++,nIter*=2)
			// if(p==0){
        	const unsigned powP = (unsigned)pow(2.0, (double)p);    

            #pragma omp parallel for
            for (int k = 0; k < (int)N / 2; ++k) {         // for butter fly calculation         // In kernel it is = thread id i.e addr i.e not N/2 rather size/2 and lThread % n makes it N/2.
                const unsigned indexAdd = i * N + (k /  powP) * 2 * powP + k % powP + offset; // getting index // offset is the location of data in real memory starting point

                const unsigned indexMult = indexAdd + powP;

                const unsigned kk = k % powP;
                
                
                if (fft_type == 0) {
                	cs = cos(twopi * kk / ( 2 * powP));          // coefficent calcultaions Wn  // at stage 0, 2*powP = 2 and kk= 0 (means W at index 0 and index 1 are same with opposite sign) and at stage =1, 2*powP = 4 and kk = 0 and 2 
                	sn = sin(twopi * kk / ( 2 * powP));          // coefficent calcultaions Wn
            	}else if (fft_type == 1) {
                	cs = cos(twopi * kk / ( 2 * powP));          // coefficent calcultaions Wn  // at stage 0, 2*powP = 2 and kk= 0 (means W at index 0 and index 1 are same with opposite sign) and at stage =1, 2*powP = 4 and kk = 0 and 2 
                	sn = -1*sin(twopi * kk / ( 2 * powP));          // coefficent calcultaions Wn
            	}


				const double addReal = data3DRr[indexAdd + ASPAN_offset];// +i%yRange+yFRange
				const double addImag = data3DRi[indexAdd + ASPAN_offset];

				// complex numbers multiplication start// W * x(1) // tempReal = W * x(1) and addReal = x(0), and indexAdd = 0 and indexMult = 1
				const double tempReal = cs * data3DRr[indexMult + ASPAN_offset] +    
				sn * data3DRi[indexMult + ASPAN_offset];
				const double tempImag = cs * data3DRi[indexMult + ASPAN_offset] - 
				sn * data3DRr[indexMult + ASPAN_offset];
				// complex numbers multiplication end

				// Complex numbers addition start
				data3DRr[indexAdd + ASPAN_offset] = addReal + tempReal;          // X(0) = x(0) + W * x(1)
				data3DRi[indexAdd + ASPAN_offset] = addImag + tempImag;
				// Complex numbers addition end

				// Complex numbers subtraction start
				data3DRr[indexMult + ASPAN_offset] =  addReal - tempReal;        // X(1) = x(0) - W * x(1)
				data3DRi[indexMult + ASPAN_offset] =  addImag - tempImag;
				// Complex numbers subtraction end

				/*	if(fft_type == 1 && p==powN-1) //if last stage
				{
					// Complex numbers addition start
					data3DRr[indexAdd + ASPAN_offset] *= norm;          // X(0) = x(0) + W * x(1)
					data3DRi[indexAdd + ASPAN_offset] *= norm;
					// Complex numbers addition end

					// Complex numbers subtraction start
					data3DRr[indexMult + ASPAN_offset] *= norm;        // X(1) = x(0) - W * x(1)
					data3DRi[indexMult + ASPAN_offset] *= norm;
					// Complex numbers subtraction end
				}*/
			}
			// }	
		}
        // for reduction
	}
    
    //const double end = omp_get_wtime();
    //totalcputime+=end - start;
   	// cout << "CPU Time: " << end - start << endl;
    /*
    if (print) {
    	printResult(size,data3DRr,data3DRi);
    }
	*/

}

/* Functions Implementation end*/
