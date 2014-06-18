#include "rw_maps.h"

#include "../fluidturtle/t_io.h"
#include "../fluidturtle/tensors3D.h"

#include "../../meteoio_plugin/meteoioplugin.h"

#include <libgen.h>
#include <sys/stat.h>
#include <sstream>

using namespace std;

/**
 * @brief copies a GeoMatrix<double> in a GeoMatrix<short>
 *
 * Performs a downcast from double to short for each value of M
 *
 * @param[out] S the matrix to fill
 * @param[in] M the matrix with the filling values
 */
void copyshort_doublematrix(GeoMatrix<short>& S, GeoMatrix<double>& M)
{

    long r, c;

    S.resize(M.getRows() , M.getCols());
    for(r = 1; r < M.getRows(); r++)
    {
        for(c = 1; c < M.getCols(); c++)
        {
            S[r][c] = (short)M[r][c];
        }
    }

}


/**
 * @brief copies a GeoMatrix<double> in a GeoMatrix<long>
 *
 * Performs a downcast from double to long for each value of M
 *
 * @param[out] L the matrix to fill
 * @param[in] M the matrix with the filling values
 */
void copylong_doublematrix(GeoMatrix<long>& L, GeoMatrix<double>& M)
{

    long r, c;

    L.resize(M.getRows(), M.getCols());
    for(r = 1; r < M.getRows(); r++)
    {
        for(c = 1; c < M.getCols(); c++)
        {
            L[r][c] = (long)M[r][c];
        }
    }

}

void copydoublematrix_const(double c0, GeoMatrix<double>& Mref, GeoMatrix<double>& M, double NOVALUE)
{

    long r, c;

    M.resize(Mref.getRows(), Mref.getCols());
    for(r = 1; r < M.getRows(); r++)
    {
        for(c = 1; c < M.getCols(); c++)
        {
            if(Mref[r][c] == NOVALUE)
            {
                M[r][c] = NOVALUE;
            }
            else
            {
                M[r][c] = c0;
            }
        }
    }

}
//----------------------

void write_suffix(std::string &suffix, long i, short start)
{
    std::stringstream lStream ;
    lStream << std::setw(4) << std::setfill('0') << i ;
    std::string lString = lStream.str();
    suffix.replace(start, lString.size(), lString);
}

std::string namefile_i(std::string name, long i)
{

    std::string SSSS = "SSSS" ;
    std::string name_out;

    write_suffix(SSSS, i, 0);

    name_out = name + SSSS + textfile ;

    return name_out;
}


std::string namefile_i_we(std::string name, long i)
{

    std::string SSSS = "LSSSS" ;
    std::string name_out;

    write_suffix(SSSS, i, 1);

    name_out = name + SSSS;

    return name_out;
}


std::string namefile_i_we2(std::string name, long i)
{

    std::string SSSS = "SSSS" ;
    std::string name_out;

    write_suffix(SSSS, i, 0);

    name_out = name + SSSS;

    return name_out ;

}

// TODO: Noori - supposed to return a pointer
GeoVector<double> read_map_vector(std::string namefile, GeoMatrix<long>& rc)
{

    GeoMatrix<double> M;
    GeoVector<double> V;
    long i, n = rc.getRows() - 1;

    meteoio_readMap(string(namefile), M);

    V.resize(n);

    for (i = 1; i < n; i++)
    {
        V[i] = M[rc[i][1]][rc[i][2]];
    }

    return V;

}



void write_map(std::string filename, short type, short format, GeoMatrix<double>& M, TInit *UV, long novalue)
{

//	type=0  floating point
//	type=1  integer

//	format=1 fluidturtle
//	format=2 grassascii
//	format=3 esriascii

    if(format == 1)
    {
        t_error("The fluidturtle format is not support any more");
    }
    else if(format == 2)
    {
        t_error("The Grass Ascii format not supported any more");
    }
    else if(format == 3)
    {
        write_esriascii(filename, type, M, UV, novalue);
    }

}

void write_map(std::string filename, short type, short format, GeoMatrix<long>& M, TInit *UV, long novalue)
{

//	type=0  floating point
//	type=1  integer

//	format=1 fluidturtle
//	format=2 grassascii
//	format=3 esriascii

    if(format == 1)
    {
        t_error("The fluidturtle format is not support any more");
    }
    else if(format == 2)
    {
        t_error("The Grass Ascii format not supported any more");
    }
    else if(format == 3)
    {
        write_esriascii(filename, type, M, UV, novalue);
    }

}

void write_map_vector(std::string filename, short type, short format, const GeoVector<double>& V, TInit *UV, long novalue, long **j, long nr, long nc)
{
    //	type=0  floating point
    //	type=1  integer

    //	format=1 fluidturtle
    //	format=2 grassascii
    //	format=3 esriascii

    if(format == 1)
    {
        t_error("The fluidturtle format is not support any more");
    }
    else if(format == 2)
    {
        t_error("Grass ascii map format not supported any more");
    }
    else if(format == 3)
    {
        write_esriascii_vector(filename, type, V, j, nr, nc, UV, novalue);
    }

}

//------------------------------------
void write_tensorseries(short a, long l, long i, std::string filename, short type, short format, GeoTensor<double>& T, TInit *UV, long novalue)
{

//	a=0 non include "l" nel suffisso
//	a=1 include "l" nel suffisso
//	l:layer
//	i:temporal step

    std::string SSSSLLLLL = "SSSSLLLLL" ;
    std::string SSSS =  "SSSS" ;
    std::string name;
    long r, c;

    GeoMatrix<double> M;

    if(a == 0)
    {
        write_suffix(SSSS, i, 0);
        name = filename ;
        name += SSSS ;
    }
    else if(a == 1)
    {
        write_suffix(SSSSLLLLL, i, 0);
        write_suffix(SSSSLLLLL, l, 5);
        name = filename;
        name += SSSSLLLLL;
    }
    else
    {
        t_error("Value not admitted");
    }


    M.resize(T.getRh(), T.getCh());
    for(r = 1; r < T.getRh(); r++)
    {
        for(c = 1; c < T.getCh(); c++)
        {
            M[r][c] = T[l][r][c];
        }
    }

    write_map(name, type, format, M, UV, novalue);

}


void write_tensorseries_vector(short a, long l, long i, std::string filename, short type, short format, GeoMatrix<double>& T, TInit *UV, long novalue, long **J, long nr, long nc)
{

    //	a=0 non include "l" nel suffisso
    //	a=1 include "l" nel suffisso
    //	l:layer
    //	i:temporal step

    std::string SSSSLLLLL = "SSSSLLLLL" ;
    std::string SSSS = "SSSS" ;

    string name;

    long j, npoints = T.getCols();

    GeoVector<double> V;

    if(a == 0)
    {
        write_suffix(SSSS, i, 0);

        name = filename + SSSS;
    }
    else if(a == 1)
    {
        write_suffix(SSSSLLLLL, i, 0);
        write_suffix(SSSSLLLLL, l, 5);

        name = filename + SSSSLLLLL;
    }
    else
    {
        t_error("Value not admitted");
    }

    V.resize(npoints + 1);

    for(j = 1; j <= npoints - 1; j++)
    {
        V[j] = T[l][j];
    }

    write_map_vector(name, type, format, V, UV, novalue, J, nr, nc);

}


//---------------------------------------------

void write_tensorseries2(std::string suf, long l, std::string filename, short type, short format, DOUBLETENSOR *T, TInit *UV, long novalue)
{

    std::string LLLLL = "LLLLL" ;
    std::string temp1, temp2;
    long r, c;

    GeoMatrix<double> M;

    temp1 = LLLLL + suf ;
    write_suffix(temp1, l, 1);

    M.resize(T->nrh + 1, T->nch + 1);

    for(r = 1; r <= T->nrh; r++)
    {
        for(c = 1; c <= T->nch; c++)
        {
            M[r][c] = T->co[l][r][c];
        }
    }

    temp2 = filename + temp1 ;
    write_map(temp2.c_str(), type, format, M, UV, novalue);

}

void write_tensorseries2_vector(std::string suf, long l, std::string filename, short type, short format, GeoMatrix<double>& T, TInit *UV, long novalue, long **J, long nr, long nc)
{

    std::string LLLLL = "LLLLL" ;
    std::string temp1, temp2;
    long i, npoints = T.getCols();
    GeoVector<double> V;

    temp1 = LLLLL + suf ;
    write_suffix(temp1, l, 1);

    V.resize(npoints + 1);

    for(i = 1; i < npoints; i++)
    {
        V[i] = T[l][i];
    }

    temp2 = filename + temp1;
    write_map_vector(temp2, type, format, V, UV, novalue, J, nr, nc);

}


//---------------
void write_tensorseries3(std::string suffix, std::string filename, short type, short format, DOUBLETENSOR *T, TInit *UV, long novalue)
{

    long l;
    for(l = T->ndl; l <= T->ndh; l++)
    {
        write_tensorseries2(suffix, l, filename, type, format, T, UV, novalue);
    }
}

//--------------------------
void write_tensorseries3_vector(std::string suffix, std::string filename, short type, short format, GeoMatrix<double>& T, TInit *UV, long novalue, long **J, long nr, long nc)
{

    long l;
    // TODO: need check nrl
    for(l = 1; l < T.getRows(); l++)
    {
        write_tensorseries2_vector(suffix, l, filename, type, format, T, UV, novalue, J, nr, nc);
    }
}


/*===================================================================================*/
/*===================functions copied from the file write_ascii.c====================*/
/*===================================================================================*/

void write_esriascii(std::string name, short type, GeoMatrix<long>& DTM, TInit *UV, long novalue)
{

//	type=0  floating point
//	type=1  integer

    FILE *f;
    long r, c;
    std::string temp;

    if(UV->U[1] != UV->U[2])
    {
        printf("\nCannot export in esriascii, grid not square, Dx=%f Dy=%f \n", UV->U[2], UV->U[1]);
        t_error("Fatal error");
    }

    temp = name + ascii_esri ;
    f = fopen(temp.c_str(), "w");

    fprintf(f, "ncols         %u\n", DTM.getCols() - 1);
    fprintf(f, "nrows         %u\n", DTM.getRows() - 1);
    fprintf(f, "xllcorner     %f\n", UV->U[4]);
    fprintf(f, "yllcorner     %f\n", UV->U[3]);
    fprintf(f, "cellsize      %f\n", UV->U[1]);
    fprintf(f, "NODATA_value  %ld.0\n", novalue);

    for(r = 1; r < DTM.getRows(); r++)
    {
        for(c = 1; c < DTM.getCols(); c++)
        {
            if((long)DTM[r][c] == novalue)
            {
                fprintf(f, "%ld.0", novalue);
            }
            else
            {
                if(type == 1)
                {
                    fprintf(f, "%ld", (long)(DTM[r][c]));
                }
                else
                {
                    fprintf(f, "%ld", DTM[r][c]);
                }
            }
            if(c < DTM.getCols()) fprintf(f, " ");
        }
        if(r < DTM.getRows()) fprintf(f, "\n");
    }
    fprintf(f, "\n"); // added by Matteo to avoid warnings when reading with R
    fclose(f);
}


void write_esriascii(std::string name, short type, GeoMatrix<double>& DTM, TInit *UV, long novalue)
{

//	type=0  floating point
//	type=1  integer

    FILE *f;
    long r, c;
    std::string temp;

    if(UV->U[1] != UV->U[2])
    {
        printf("\nCannot export in esriascii, grid not square, Dx=%f Dy=%f \n", UV->U[2], UV->U[1]);
        t_error("Fatal error");
    }

    temp = name ;
    temp += ascii_esri;
    f = fopen(temp.c_str(), "w");

    fprintf(f, "ncols         %u\n", DTM.getCols() - 1);
    fprintf(f, "nrows         %u\n", DTM.getRows() - 1);
    fprintf(f, "xllcorner     %f\n", UV->U[4]);
    fprintf(f, "yllcorner     %f\n", UV->U[3]);
    fprintf(f, "cellsize      %f\n", UV->U[1]);
    fprintf(f, "NODATA_value  %ld.0\n", novalue);

    for(r = 1; r < DTM.getRows(); r++)
    {
        for(c = 1; c < DTM.getCols(); c++)
        {
            if((long)DTM[r][c] == novalue)
            {
                fprintf(f, "%ld.000", novalue);
            }
            else
            {
                if(type == 1)
                {
                    fprintf(f, "%ld", (long)(DTM[r][c]));
                }
                else
                {
                    fprintf(f, "%.3f", DTM[r][c]);
                }
            }
            if(c < DTM.getCols()) fprintf(f, " ");
        }
        if(r < DTM.getRows()) fprintf(f, "\n");
    }
    fprintf(f, "\n"); // added by Matteo to avoid warnings when reading with R
    fclose(f);
}

//===============

void write_esriascii_vector(char *name, short type, const GeoVector<double>& DTM, long **j, long nr, long nc, TInit *UV, long novalue)
{
    //	type=0  floating point
    //	type=1  integer

    FILE *f;
    long r, c;
    std::string temp;

    if(UV->U[1] != UV->U[2])
    {
        printf("\nCannot export in esriascii, grid not square, Dx=%f Dy=%f \n", UV->U[2], UV->U[1]);
        t_error("Fatal error");
    }

    temp = name ;
    temp += ascii_esri;
    f = fopen(temp.c_str(), "w");

    fprintf(f, "ncols         %ld\n", nc);
    fprintf(f, "nrows         %ld\n", nr);
    fprintf(f, "xllcorner     %f\n", UV->U[4]);
    fprintf(f, "yllcorner     %f\n", UV->U[3]);
    fprintf(f, "cellsize      %f\n", UV->U[1]);
    fprintf(f, "NODATA_value  %ld.0\n", novalue);

    for(r = 1; r <= nr; r++)
    {
        for(c = 1; c <= nc; c++)
        {
            if (j[r][c] > 0)
            {
                if(type == 1)
                {
                    fprintf(f, "%ld", (long)(DTM[j[r][c]]));
                }
                else
                {
                    fprintf(f, "%.3f", DTM[j[r][c]]);
                }
            }
            else
            {
                fprintf(f, "%ld.000", novalue);
            }
            if(c < nc) fprintf(f, " ");
        }
        if(r < nr) fprintf(f, "\n");
    }
    fclose(f);
}

void write_esriascii_vector(string name, short type, const GeoVector<double>& DTM, long **j, long nr, long nc, TInit *UV, long novalue)
{
    //	type=0  floating point
    //	type=1  integer

    char *basedir;
    int ret = 0;

    FILE *f;
    long r, c;
    string temp;

    if(UV->U[1] != UV->U[2])
    {
        printf("\nCannot export in esriascii, grid not square, Dx=%f Dy=%f \n", UV->U[2], UV->U[1]);
        t_error("Fatal error");
    }

    temp = name + ascii_esri;

    char * lStrBase = strdup(temp.c_str()) ;
    basedir = dirname(lStrBase);
    ret = mkdirp(basedir, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    free(lStrBase) ;
    if(-1 == ret)
    {
        t_error("write_esriascii_vector(): Unable to create parent directories of file" + temp);
    }
    f = fopen(temp.c_str(), "w");
    if(NULL == f)
    {
        t_error("write_esriascii_vector(): Unable to open file `" + temp + "` in write mode.");
    }

    fprintf(f, "ncols         %ld\n", nc);
    fprintf(f, "nrows         %ld\n", nr);
    fprintf(f, "xllcorner     %f\n", UV->U[4]);
    fprintf(f, "yllcorner     %f\n", UV->U[3]);
    fprintf(f, "cellsize      %f\n", UV->U[1]);
    fprintf(f, "NODATA_value  %ld.0\n", novalue);

    for(r = 1; r <= nr; r++)
    {
        for(c = 1; c <= nc; c++)
        {
            if (j[r][c] > 0)
            {
                if(type == 1)
                {
                    fprintf(f, "%ld", (long)(DTM[j[r][c]]));
                }
                else
                {
                    fprintf(f, "%.3f", DTM[j[r][c]]);
                }
            }
            else
            {

                fprintf(f, "%ld.000", novalue);
            }
            if(c < nc) fprintf(f, " ");
        }
        if(r < nr) fprintf(f, "\n");
    }
    fclose(f);


}


//-----------------------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------------------

void error_message(short format, long n, long n1, long n2, long n3, char *name)
//format=1 grassascii
//format=2 esriascii

{
    if(n == n1 || n == n2 || n == n3)
    {
        std::string lFile0 = name ;
        lFile0 += ascii_grass ;

        std::string lFile1 = name ;
        lFile1 += ascii_esri ;

        if(format == 1) printf("File -3- %s incompleted, end of file or end of line reached", lFile0.c_str());
        if(format == 2) printf("File -4- %s incompleted, end of file or end of line reached", lFile1.c_str());
        t_error("Fatal error");
    }
}


/*===================functions copied from geomorphology.0875.c============*/

void curvature(double deltax, double deltay, GeoMatrix<double>& topo, GeoMatrix<double>& c1, GeoMatrix<double>& c2, GeoMatrix<double>& c3, GeoMatrix<double>& c4, long undef)
{

    long r, c;
    long R1, R2, C1, C2;
    long nc = topo.getCols() - 1;
    long nr = topo.getRows() - 1;
    double delta;

    // Compute the curvature.
    for(r = 1; r <= nr; r++)
    {
        for(c = 1; c <= nc; c++)
        {
            if((long)topo[r][c] != undef)
            {

                c1[r][c] = 0.0;
                R1 = r - 1;
                R2 = r + 1;
                C1 = c;
                C2 = c;
                delta = deltay;
                if(R1 >= 1 && R1 <= nr && R2 >= 1 && R2 <= nr && C1 >= 1 && C1 <= nc && C2 >= 1 && C2 <= nc)
                {
                    if((long)topo[R1][C1] != undef && (long)topo[R2][C2] != undef)
                    {
                        c1[r][c] += (topo[R1][C1] + topo[R2][C2] - 2.*topo[r][c]) / pow(delta, 2.);
                    }
                }

                c2[r][c] = 0.0;
                R1 = r;
                R2 = r;
                C1 = c + 1;
                C2 = c - 1;
                delta = deltax;
                if(R1 >= 1 && R1 <= nr && R2 >= 1 && R2 <= nr && C1 >= 1 && C1 <= nc && C2 >= 1 && C2 <= nc)
                {
                    if((long)topo[R1][C1] != undef && (long)topo[R2][C2] != undef)
                    {
                        c2[r][c] += (topo[R1][C1] + topo[R2][C2] - 2.*topo[r][c]) / pow(delta, 2.);
                    }
                }

                c3[r][c] = 0.0;
                R1 = r - 1;
                R2 = r + 1;
                C1 = c - 1;
                C2 = c + 1;
                delta = sqrt(deltax * deltay);
                if(R1 >= 1 && R1 <= nr && R2 >= 1 && R2 <= nr && C1 >= 1 && C1 <= nc && C2 >= 1 && C2 <= nc)
                {
                    if((long)topo[R1][C1] != undef && (long)topo[R2][C2] != undef)
                    {
                        c3[r][c] += (topo[R1][C1] + topo[R2][C2] - 2.*topo[r][c]) / pow(delta, 2.);
                    }
                }

                c4[r][c] = 0.0;
                R1 = r - 1;
                R2 = r + 1;
                C1 = c + 1;
                C2 = c - 1;
                delta = sqrt(deltax * deltay);
                if(R1 >= 1 && R1 <= nr && R2 >= 1 && R2 <= nr && C1 >= 1 && C1 <= nc && C2 >= 1 && C2 <= nc)
                {
                    if((long)topo[R1][C1] != undef && (long)topo[R2][C2] != undef)
                    {
                        c4[r][c] += (topo[R1][C1] + topo[R2][C2] - 2.*topo[r][c]) / pow(delta, 2.);
                    }
                }


            }
            else
            {

                c1[r][c] = (double)undef;
                c2[r][c] = (double)undef;
                c3[r][c] = (double)undef;
                c4[r][c] = (double)undef;

            }

        }
    }
}

//------------------------
short is_boundary(long r, long c, GeoMatrix<double>& dem, long novalue)
{

    long ir, ic;
    short yes = 0;

    ir = -1;
    ic = 0;
    if( (long)dem[r + ir][c + ic] == novalue ) yes = 1;

    ir = -1;
    ic = 1;
    if( (long)dem[r + ir][c + ic] == novalue ) yes = 1;

    ir = 0;
    ic = 1;
    if( (long)dem[r + ir][c + ic] == novalue ) yes = 1;

    ir = 1;
    ic = 1;
    if( (long)dem[r + ir][c + ic] == novalue ) yes = 1;

    ir = 1;
    ic = 0;
    if( (long)dem[r + ir][c + ic] == novalue ) yes = 1;

    ir = 1;
    ic = -1;
    if( (long)dem[r + ir][c + ic] == novalue ) yes = 1;

    ir = 0;
    ic = -1;
    if( (long)dem[r + ir][c + ic] == novalue ) yes = 1;

    ir = -1;
    ic = -1;
    if( (long)dem[r + ir][c + ic] == novalue ) yes = 1;

    return yes;

}


long row(double N, long nrows, TInit *UV, long novalue)
{

    long cnt;

    if(N < UV->U[3] || N > UV->U[3] + nrows * UV->U[1])
    {
        return novalue;
    }
    else
    {
        cnt = 0;
        do
        {
            cnt++;
        }
        while(UV->U[3] + (nrows - cnt)*UV->U[1] > N);
        return cnt;
    }
}


long col(double E, long ncols, TInit *UV, long novalue)
{

    long cnt;

    if(E < UV->U[4] || E > UV->U[4] + ncols * UV->U[2])
    {
        return novalue;
    }
    else
    {
        cnt = 0;
        do
        {
            cnt++;
        }
        while(UV->U[4] + cnt * UV->U[2] < E);
        return cnt;
    }
}



//Presa da geomorphology099 e modificato der_min
void nablaquadro_mask(GeoMatrix<double>& Z0, GeoMatrix<short>& curv, GeoVector<double>& U, GeoVector<double>& V)

{

    short y;
    long i, j, h, rows, cols;
    double grid[9], z[9], derivate2;
    double der_min = 0.00001; /*limite per la limite per la planarita'*/

    short v[13][2] = {
        { 0, 0},
        { 0, 1},
        { -1, 1},
        { -1, 0},
        { -1, -1},
        { 0, -1},
        { 1, -1},
        { 1, 0},
        { 1, 1},
        { 0, 0},
        { 0, 0},
        { 0, 0},
        { 0, 0}
    };

    grid[0] = 0;
    grid[1] = grid[5] = U[1];
    grid[3] = grid[7] = U[2];
    grid[2] = grid[4] = grid[6] = grid[8] = sqrt(grid[1] * grid[1] + grid[3] * grid[3]);

    rows = Z0.getRows() - 1;
    cols = Z0.getCols() - 1;

    for(i = 2; i <= rows - 1; i++)
    {
        for(j = 2; j <= cols - 1; j++)
        {
            z[0] = Z0[i][j];
            if(z[0] != V[2])
            {
                y = 1;
                for(h = 1; h <= 8; h++)
                {
                    z[h] = Z0[i + v[h][0]][j + v[h][1]];
                    if(z[h] == V[2])
                    {
                        y = 0;
                        break;
                    }
                }
                if(y == 0)
                {
                    curv[i][j] = 1;
                }
                else
                {
                    derivate2 = 0.5 * ((z[1] + z[5] - 2 * z[0]) / (grid[1] * grid[1]) + (z[3] + z[7] - 2 * z[0]) / (grid[3] * grid[3]));
                    derivate2 = derivate2 + 0.5 * ((z[2] + z[4] + z[6] + z[8] - 4 * z[0]) / (grid[6] * grid[6]));

                    if(fabs(derivate2) <= der_min || derivate2 > der_min) //plane or concave
                    {
                        curv[i][j] = 0;
                    }
                    else
                    {
                        curv[i][j] = 1;		//convex
                    }
                }
            }
        }
    }
}

/*====================copied function from init.c ==================*/

void initmatrix(double val, GeoMatrix<double>& destination, GeoMatrix<double>& origin, double novalue)
{

    long r, c;
    for(r = 1; r < destination.getRows(); r++)
    {
        for(c = 1; c < destination.getCols(); c++)
        {
            if(origin[r][c] != novalue) destination[r][c] = val;
        }
    }
}



//----------------------------------
