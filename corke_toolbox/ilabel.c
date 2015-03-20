/* ilabel.c
 *
 * Fast image labelling funciton.
 *
 *	L = ILABEL(IM [,OPT])
 *	[L,LMAX] = ILABEL(IM)
 *	[L,LMAX,PARENTS] = ILABEL(IM)
 *
 *	$Header: /home/autom/pic/cvsroot/image-toolbox/ilabel.c,v 1.5 2005/11/24 11:14:36 pic Exp $
 *
 *	$Log: ilabel.c,v $
 *	Revision 1.5  2005/11/24 11:14:36  pic
 *	Change to copyright notices.
 *
 *	Revision 1.4  2005/10/30 02:52:33  pic
 *	Use old style comments, tidyup error message.
 *
 *	Revision 1.3  2005/07/03 10:47:52  pic
 *	Add support for 4 and 8-way connectivity.
 *
 *	Revision 1.2  2004/12/03 07:30:38  pic
 *	Removed debug.
 *	
 *	Revision 1.1  2002/08/28 04:53:06  pic
 *	Initial CVS version.
 *	
 *	Revision 1.2  2000/03/10 07:04:11  pic
 *	Added region hierarchy and return of region parent vector.
 *
 *
 *   By Peter Corke, Copyright (c) CSIRO, 1995  Machine Vision Toolbox for Matlab
 *		pic 3/95
 */
#include "mex.h"
#include <math.h>

/*
#define	DEBUG
*/

/* Input Arguments */

#define	IM_IN		prhs[0]
#define	OPT_IN		prhs[1]

/* Output Arguments */

#define	IM_OUT	plhs[0]
#define	MAX_OUT	plhs[1]
#define	PARENT_OUT	plhs[2]

int	nrows, ncols;
#define	MAXLABEL	10000
#define	UNKNOWN		-1

#define	PIX(v,r,c)	v[(r)+(c)*nrows]
int lresolve(int label);

static int	*lmap;
int		connectivityWays;

/*
 * actual labelling code:
 *
 *	im	double prec. input image (input)
 *	dlim	double prec. label image (output)
 */
static int 
ilabel(double *im, double *dlim, int **Parents)
{
	int	*lim, *lmap2, row, col, i, j, k, nlabels;
	int	prevlab, curlab, curpix, prevpix;
	int	newlabel;
	int	*parents;

	/* allocate label map and initialize to zero */
	lmap = mxCalloc(MAXLABEL, sizeof(int));
	lmap2 = mxCalloc(MAXLABEL, sizeof(int));

	parents = mxCalloc(MAXLABEL, sizeof(int));

	/* tempory integer labelled image (for speed) */
	lim = mxCalloc(nrows*ncols, sizeof(int));


	/*
	 * first pass labelling loop.  Only does 4-way connectivity
	 */
	newlabel = 0;
	for (row=0; row<nrows; row++) {
		prevlab = UNKNOWN;
		for (col=0; col<ncols; col++) {
			curpix = PIX(im,row,col);

			/* if no change in pixel value then inherit label from left */
			curlab = UNKNOWN;
			if ((col > 0) && (curpix == prevpix))
				curlab = prevlab;
			
			/*
			 * check whether a label merge should happen, adjacent
			 * pixels with the same value but different labels
			 * means that one should change.
			 *
			 * merge can only happen on second row onwards 
			 */
			if (row > 0) {
				if (	(PIX(im,row-1,col) == curpix) &&
					(lresolve(PIX(lim,row-1,col)) != curlab)
				) {
					/* we have a merge to N */
					int	newlabel;
					
					newlabel = lresolve(PIX(lim,row-1,col));
					/*
					newlabel = PIX(lim,row-1,col);
					*/

#ifdef	DEBUG
					printf("mergeN(%d,%d): %d becomes %d\n",
						row, col, curlab, newlabel);
#endif
					if (curlab != UNKNOWN)
						lmap[curlab] = newlabel;
					curlab = newlabel;

				} else if (
					connectivityWays == 8 &&
					(col > 0) &&
					(PIX(im,row-1,col-1) == curpix) &&
					(lresolve(PIX(lim,row-1,col-1)) != curlab)
				) {
					/* we have a merge to NW */
					int	newlabel;
					
					newlabel = lresolve(PIX(lim,row-1,col-1));
					/*
					newlabel = PIX(lim,row-1,col);
					*/

#ifdef	DEBUG
					printf("mergeNW(%d,%d): %d becomes %d\n",
						row, col, curlab, newlabel);
#endif
					if (curlab != UNKNOWN)
						lmap[curlab] = newlabel;
					curlab = newlabel;

				} else if (
					connectivityWays == 8 &&
					(col < (ncols-1)) &&
					(PIX(im,row-1,col+1) == curpix) &&
					(lresolve(PIX(lim,row-1,col+1)) != curlab)
				) {
					/* we have a merge to NE */
					int	newlabel;
					
					newlabel = lresolve(PIX(lim,row-1,col+1));
					/*
					newlabel = PIX(lim,row-1,col);
					*/

#ifdef	DEBUG
					printf("mergeNE(%d,%d): %d becomes %d\n",
						row, col, curlab, newlabel);
#endif
					if (curlab != UNKNOWN)
						lmap[curlab] = newlabel;
					curlab = newlabel;

				}

			}

			if ((row > 0) && (col > 0)) {
				/*
				 * check for enclosure
				 */
				int	left, above, diag;

				left = prevlab;
				above = lresolve( PIX(lim,row-1,col) );
				diag = lresolve( PIX(lim,row-1,col-1) );
				if (	(left == curlab) &&
					(above == curlab) &&
					(diag != curlab)
				) {
#ifdef	DEBUG
					printf("label %d encloses %d\n",
						curlab, i);
#endif
					/* we have an enclosure */
					parents[diag] = curlab;
				}
			}

			/* if label still not known, assign new */
			if (curlab == UNKNOWN) {
				curlab = ++newlabel;
				if (newlabel >= MAXLABEL)
					mexErrMsgTxt("ilabel: too many regions");
#ifdef	DEBUG
				printf("new label(%d,%d): %d\n", 
					row, col, curlab);
#endif
			}

			PIX(lim,row,col) = curlab;
			prevlab = curlab;
			prevpix = curpix;
		}
	}

#ifdef	DEBUG
	printf("max lim is %d\n", newlabel);
#endif

	/*
	 * now eliminate redirections from the label map
	 *
	 *	lmap[pass1 label] -> pass 2 label
	 *	lmap2[pass 2 label] -> final label
	 */
#ifdef	DEBUG
	printf("----------------------\nlmap:\n");
#endif
	for (i=1,nlabels=0; i<=newlabel; i++) {
#ifdef	DEBUG
		printf("(%d) = %d\n", i, lmap[i]);
#endif
		if (lmap[i] == 0)
			lmap2[i] = ++nlabels;	/* assign new sequential label */
	}

	/*
	 * now adjust the label map so that consecutive labels appear in the
	 * labelled image, ie. no missing labels.
	 */
	for (i=0; i<=newlabel; i++)
		if (lmap[i] != 0) {
			j = lresolve(i);
			lmap2[i] = lmap2[j];
		}
#ifdef	DEBUG
	printf("----------------------\nlmap2:\n");
	for (i=1; i<=newlabel; i++)
		printf("(%d) = %d\n", i, lmap2[i]);
#endif

#ifdef	DEBUG
	printf("----------------------\nparents:\n");
	for (i=1; i<=newlabel; i++)
		printf("parent[%d] = %d\n", i, parents[i]);
#endif

	/*
	 * resolve the labels in the parent array and assign to double proc
	 * output array
	 */
	*Parents = mxCalloc(nlabels, sizeof(int));
	for (i=0; i<=newlabel; i++) {
		int	par = parents[i], child;

		if (par) {
			child = lmap2[i];
			par = lmap2[par];
			(*Parents)[child-1] = par;
		}
	}
	mxFree(parents);

	/*
	 * resolve the labels in the integer labelled image and assign
	 * to the double prec. output image
	 */
	for (row=0; row<nrows; row++)
		for (col=0; col<ncols; col++)
			PIX(dlim,row,col) = lmap2[ PIX(lim,row,col) ];

	mxFree(lmap);
	mxFree(lmap2);
	mxFree(lim);

	return(nlabels);
}

/*
 * resolve a label to it's true value via the label map
 */
int
lresolve(int l)
{
	int	i;

	for (i=l; lmap[i] > 0; )
		i = lmap[i];
#ifdef	DEBUG
	printf("resolved %d to %d\n", l, i);
#endif
	return i;
}


/*
 * MATLAB interface function, check arguments, then call ilabel() above
 */
void
mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	mxArray	*Mlabel, *Mmax;
	double	*im, *label;
	int	maxlabel, *parents;

	/* Check for proper number of arguments */

	connectivityWays = 4;
	switch (nrhs) {
	case 2: {
		double	opt = mxGetScalar(OPT_IN);

		if (opt == 8)
			connectivityWays = 8;
		else if (opt != 4)
			mexErrMsgTxt("ILABEL connectivity must be 4 or 8.");
		}
	case 1:
		nrows = mxGetM(IM_IN);
		ncols = mxGetN(IM_IN);
		/*
		printf("size %d x %d\n", nrows, ncols);
		*/
		if (!mxIsNumeric(IM_IN) || mxIsComplex(IM_IN) || 
			!mxIsDouble(IM_IN)) {
			mexErrMsgTxt("ILABEL requires a real matrix.");
		}
		break;
	default:
		mexErrMsgTxt("ILABEL requires one or more input arguments");
		break;
	}



	/* Create a matrix for the return argument */
	Mlabel = mxCreateDoubleMatrix(nrows, ncols, mxREAL);

	im = mxGetPr(IM_IN);
	label = mxGetPr(Mlabel);

	/* Do the actual computations in a subroutine */
	maxlabel = ilabel(im, label, &parents);


	switch (nlhs) {
	case 3: {
		double	*p;
		int	i;

		PARENT_OUT = mxCreateDoubleMatrix(maxlabel, 1, mxREAL);
		p = mxGetPr(PARENT_OUT);
		for (i=0; i<maxlabel; i++)
			p[i] = parents[i];
	    }
		/* fall through */
	case 2: {
		double	*p;

		Mmax = mxCreateDoubleMatrix(1, 1, mxREAL);
		p = mxGetPr(Mmax);
		*p = maxlabel;
		MAX_OUT = Mmax;
		}
		/* fall through */
	case 1:
		IM_OUT = Mlabel;
		break;
	}


	return;
}
