#ifndef FILE_OUTPUT_HH
#define FILE_OUTPUT_HH

#include <cstdio>
#include <cstdlib>
/** Writes binary data to a file. If an error was encountered, the routine
 * prints an error message and exits.
 * \param[in] bu a pointer to the binary data to write.
 * \param[in] size the size in bytes of each record in the data.
 * \param[in] count the number records to write.
 * \param[in] fp the file handle to write to. */
inline void sfwrite(const void *bu,size_t size,size_t count,FILE *fp) {
	if(fwrite(bu,size,count,fp)!=count) {
		fputs("Error writing data to file",stderr);
		exit(1);
	}
}

/** Outputs a two dimensional array in the Gnuplot binary format. The routine
 * treats the array as periodic and repeats the first row and column.
 * \param[in] fp a file handle to write to.
 * \param[in] fld the array to use.
 * \param[in] m the horizontal grid size.
 * \param[in] n the vertical grid size.
 * \param[in] (ax,bx) the horizontal range to use.
 * \param[in] (ay,by) the vertical range to use. */
template<class T>
void gnuplot_output(FILE *fp,T *fld,int m,int n,double dx,double dy,double ax,double bx,double ay,double by) {
	int i,j;

	// Write header line
	float *fbuf=new float[m+1],*ff=fbuf;
	T *pp;
	*(ff++)=static_cast<float>(m);for(i=0;i<m;i++) *(ff++)=static_cast<float> (ax+i*dx);
	sfwrite(fbuf,sizeof(float),m+1,fp);

	// Write field entries line-by-line
	for(j=0;j<n;j++) {

		// Write header entry
		ff=fbuf;*(ff++)=static_cast<float>(ay+j*dy);

		// Write a horizontal line to the buffer
		pp=fld+j*m;for(i=0;i<m;i++) *(ff++)=static_cast<float>(*(pp++));
		sfwrite(fbuf,sizeof(float),m+1,fp);
	}

	// Remove temporary memory and close file
	delete [] fbuf;
}
#endif
