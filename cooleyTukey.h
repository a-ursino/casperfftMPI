#ifndef __COOLEYTUKEY__
#define __COOLEYTUKEY__


// Functions Prototype 

bool runCooleyTukey(const char* const argv[] );
void transpose2(double *data3DRr,double *data3DRi, int xR, int yR, int off);
void transpose3(const unsigned offset, double *data3DRr, double *data3DRi, int zR, int yR, int xR, int off);
void convolveCPU(const unsigned offset1,int ASPAN);
void cooleyTukeyCpu3DFFT(const unsigned offset, const unsigned  N, const unsigned size,double *data3DFr,double *data3DFi,double *data3DRr,double *data3DRi,int ASPAN_offset, int show_results, int fft_type);
void cooleyTukeyCpu(const unsigned offset, const unsigned  N, const unsigned size,double *data3DFr,double *data3DFi,double *data3DRr,double *data3DRi,int ASPAN_offset, int show_results, int fft_type, int full, int red, int reduction);

// Functions Prototype end

#endif
