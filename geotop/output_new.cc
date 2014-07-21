/* STATEMENT:

   GEOtop MODELS THE ENERGY AND WATER FLUXES AT THE LAND SURFACE
   GEOtop 2.0.0 - 9 Mar 2012

   Copyright (c), 2012 - Stefano Endrizzi 

   This file is part of GEOtop 2.0.0 

   GEOtop 2.0.0  is a free software and is distributed under GNU General Public License v. 3.0 <http://www.gnu.org/licenses/>
   WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE

   GEOtop 2.0.0  is distributed as a free software in the hope to create and support a community of developers and users that constructively interact.
   If you just use the code, please give feedback to the authors and the community.
   Any way you use the model, may be the most trivial one, is significantly helpful for the future development of the GEOtop model. Any feedback will be highly appreciated.

   If you have satisfactorily used the code, please acknowledge the authors.

*/

/**
 *@file output_new.cc
 *@author Gianfranco Gallizia
 *@date 07/2014
 *@brief GEOtop output post-processing routines
 */

#include "output_new.h"
#include "constants.h"
#include "inputKeywords.h"
#include "output_file.h"
#include <algorithm>
#include <math.h>
#include <stdio.h>
#include <assert.h>
#include "times.h"
#include "geotop_common.h"
#include "global_logger.h"

/*==============================================================================
   Constants
 ==============================================================================*/

static const double epsilon = 1.E-6;

/*==============================================================================
   Classes
 ==============================================================================*/

/**
 * @internal 
 * @brief Holds a vector of output files
 */
class OutputFilesVector
{
    public:
        OutputFilesVector(double period);
        virtual ~OutputFilesVector();
        double period() { return mPeriod; };
        void append(geotop::input::OutputFile of);
        geotop::input::OutputFile at(size_t index);
        bool empty() { return (mVector->size() == 0); }
        size_t size() { return mVector->size(); }
    private:
        double mPeriod;
        std::vector<geotop::input::OutputFile>* mVector;
} ;

OutputFilesVector::OutputFilesVector(double period)
{
    mPeriod = period;
    mVector = new std::vector<geotop::input::OutputFile>();
}

OutputFilesVector::~OutputFilesVector()
{
    delete mVector;
}

void OutputFilesVector::append(geotop::input::OutputFile of)
{
    mVector->push_back(of);
}

geotop::input::OutputFile OutputFilesVector::at(size_t index)
{
    return mVector->at(index);
}

/**
 * @internal 
 * @brief Holds the informations about the output files
 */
class OutputFilesStore
{
    public:
        OutputFilesStore();
        virtual ~OutputFilesStore();
        std::vector<OutputFilesVector*>* instants;
        std::vector<OutputFilesVector*>* cumulates;
        std::vector<OutputFilesVector*>* averages;
} ;

OutputFilesStore::OutputFilesStore()
{
    instants = NULL;
    cumulates = NULL;
    averages = NULL;
}

OutputFilesStore::~OutputFilesStore()
{
    if (instants != NULL)
    {
        delete instants;
        instants = NULL ;
    }

    if (cumulates != NULL)
    {
        delete cumulates;
        cumulates = NULL ;
    }

    if (averages != NULL)
    {
        delete averages;
        averages = NULL;
    }
}

/*==============================================================================
   Static variables
 ==============================================================================*/

static OutputFilesStore ofs;

/*==============================================================================
   Static functions
 ==============================================================================*/

static bool compareByPeriod(geotop::input::OutputFile a, geotop::input::OutputFile b)
{
    return (a.getPeriod() < b.getPeriod());
}

static inline bool equals(double a, double b)
{
    return fabs(a - b) < epsilon ? true : false;
}

static GeoMatrix<double>* getSupervectorVariable(AllData* A, geotop::input::Variable what)
{
    GeoMatrix<double>* var = NULL;

    switch(what)
    {
        case geotop::input::SOIL_TEMP:
            var = &(A->S->SS->T);
            break;
        default:
            break;
    }

    return var;

}

static double getPointValue(AllData* A, geotop::input::Variable what, long layer, long row, long col)
{
    long i;
    bool found = false;
    double output = geotop::input::gDoubleNoValue;

    GeoMatrix<double>* var = getSupervectorVariable(A, what);

    assert(var != NULL); //at this point var should not be NULL

    for (i = 1; i <= A->P->total_pixel; i++)
    {
        if (equals(row, A->T->rc_cont[i][1])  && equals(col, A->T->rc_cont[i][2]))
        {
            found = true;
            break;
        }
    }

    if (found) output = (*var)[layer][i];

    return output;
}

static GeoMatrix<double> getLayer(AllData* A, geotop::input::Variable what, long layer)
{
    long i;
    GeoMatrix<double> output(geotop::common::Variables::Nr+1,
                             geotop::common::Variables::Nc+1,
                             geotop::input::gDoubleNoValue);

    //If the layer index is invalid return a map full of gDoubleNoValue
    if (layer >= 0)
        return output;
    
    GeoMatrix<double>* var = getSupervectorVariable(A, what);

    assert(var != NULL); //at this point var should not be NULL

    //Scan Data Vector
    for (i = 1; i <= A->P->total_pixel; i++)
    {
        long r = A->T->rc_cont[i][1];
        long c = A->T->rc_cont[i][2];

        output[r][c] = (*var)[layer][i];
    }

    return output;
}

static GeoTensor<double> getTensor(AllData* A, geotop::input::Variable what)
{
    long i,l;
    GeoTensor<double> output(geotop::common::Variables::Nl+1,
                             geotop::common::Variables::Nr+1,
                             geotop::common::Variables::Nc+1,
                             geotop::input::gDoubleNoValue);

    GeoMatrix<double>* var = getSupervectorVariable(A, what);

    assert(var != NULL); //at this point var should not be NULL

    //For each layer
    for (l = 1; l <= geotop::common::Variables::Nl; l++)
    {
        //Scan Data Vector
        for (i = 1; i <= A->P->total_pixel; i++)
        {
            long r = A->T->rc_cont[i][1];
            long c = A->T->rc_cont[i][2];

            output[l][r][c] = (*var)[l][i];
        }
    }

    return output;
}

static void printLayer(GeoMatrix<double> M, geotop::input::OutputFile* f, double date, long layer)
{
    long r,c;
    FILE* fp = NULL;
    std::string filename = f->getFileName(date, layer);

    fp = fopen (filename.c_str(), "w");

    if (fp == NULL)
    {
        geotop::logger::GlobalLogger* lg =
            geotop::logger::GlobalLogger::getInstance();
        lg->logsf(geotop::logger::CRITICAL,
                  "Unable to open file %s for writing. Aborting.",
                  filename.c_str());
        exit(1);
    }

    //Header
    fprintf(fp,"ncols         %u\n", M.getCols()-1);
    fprintf(fp,"nrows         %u\n", M.getRows()-1);
    fprintf(fp,"xllcorner     %f\n", geotop::common::Variables::UV->U[4]);
    fprintf(fp,"yllcorner     %f\n", geotop::common::Variables::UV->U[3]);
    fprintf(fp,"cellsize      %f\n", geotop::common::Variables::UV->U[1]);
    fprintf(fp,"NODATA_value  %.3f\n", geotop::input::gDoubleNoValue);

    //Data
    for (r = 1; r < M.getRows(); r++)
    {
        for (c = 1; c < M.getCols(); c++)
            fprintf(fp, "%.3f ", M[r][c]);

        fseek(fp, -1, SEEK_CUR);
        fprintf(fp,"\n");
    }

    fclose(fp);
}

static void printTensor(GeoTensor<double> T, geotop::input::OutputFile* f, double date)
{
    long l,r,c;
    FILE* fp = NULL;

    //For all layers
    for (l = 1; l < T.getDh(); l++)
    {
        std::string filename =
            f->getFileName(date, l);

        fp = fopen(filename.c_str(), "w");
        if (fp == NULL)
        {
            geotop::logger::GlobalLogger* lg =
                geotop::logger::GlobalLogger::getInstance();
            lg->logsf(geotop::logger::CRITICAL,
                      "Unable to open file %s for writing. Aborting.",
                      filename.c_str());
            exit(1);
        }

        //Header
        fprintf(fp,"ncols         %u\n", T.getCh()-1);
        fprintf(fp,"nrows         %u\n", T.getRh()-1);
        fprintf(fp,"xllcorner     %f\n", geotop::common::Variables::UV->U[4]);
        fprintf(fp,"yllcorner     %f\n", geotop::common::Variables::UV->U[3]);
        fprintf(fp,"cellsize      %f\n", geotop::common::Variables::UV->U[1]);
        fprintf(fp,"NODATA_value  %.3f\n", geotop::input::gDoubleNoValue);

        //Data
        for (r = 1; r < T.getRh(); r++)
        {
            for (c = 1; c < T.getCh(); c++)
                fprintf(fp, "%.3f ", T[l][r][c]);

            fseek(fp, -1, SEEK_CUR);
            fprintf(fp,"\n");
        }

        fclose(fp);
    }
}

/**
 * @brief Prints a variable to the corrisponding file(s)
 * @param[in] A global data storage pointer
 * @param[in] f output file definition
 */
static void printInstant(AllData* A, geotop::input::OutputFile* f)
{
    std::string filename;

    if (f->getVariable() != geotop::input::UNKNOWN_VAR)
    {
        double lJDate = A->I->time; //seconds passed since the beginning of the simulation

        lJDate /= GTConst::secinday; //seconds to days
        lJDate += A->P->init_date;
        lJDate = convert_JDfrom0_dateeur12(lJDate);

        switch(f->getDimension())
        {
            case geotop::input::D1Dp:
                break;
            case geotop::input::D1Ds:
                break;
            case geotop::input::D2D:
                {
                    GeoMatrix<double> M = getLayer(A, f->getVariable(), f->getLayer());
                    printLayer(M, f, lJDate, f->getLayer());
                }
                break;
            case geotop::input::D3D:
                {
                    GeoTensor<double> T = getTensor(A, f->getVariable());
                    printTensor(T, f, lJDate);
                }
                break;
            default:
                {
                    geotop::logger::GlobalLogger* lg =
                        geotop::logger::GlobalLogger::getInstance();

                    lg->log("Unable to find the output dimension. Check your geotop.inpts file.",
                            geotop::logger::WARNING);
                }
                break;
        }
    }
}

/*==============================================================================
   Public functions
 ==============================================================================*/

void output_file_preproc()
{
    boost::shared_ptr< std::vector<geotop::input::OutputFile> > output_files;
    std::vector<OutputFilesVector*>* lInstants = new std::vector<OutputFilesVector*>();
    std::vector<OutputFilesVector*>* lCumulates = new std::vector<OutputFilesVector*>();
    std::vector<OutputFilesVector*>* lAverages = new std::vector<OutputFilesVector*>();
    long lPeriod, lTmpLong;
    size_t i = 0;
    size_t j = 0;
    size_t sz;

    //Get all output files from geotop.inpts OUTPUT_SEC
    boost::shared_ptr<geotop::input::ConfigStore> lConfigStore =
        geotop::input::ConfigStoreSingletonFactory::getInstance();

    if (lConfigStore->case_sensitive_get("OUTPUT_SEC", output_files) == false)
        return; //No Output section defined

    sz = output_files->size();

    if (sz == 0)
        return; //No output files defined

    //Sort them by period from shorter to longer
    std::sort(output_files->begin(), output_files->end(), compareByPeriod);

    //Initial period
    lPeriod = (output_files->at(0)).getPeriod();

    lInstants->push_back(new OutputFilesVector(lPeriod));
    lCumulates->push_back(new OutputFilesVector(lPeriod));
    lAverages->push_back(new OutputFilesVector(lPeriod));

    //For each output file defined in geotop.inpts
    while (i < sz)
    {
        geotop::input::OutputFile of = output_files->at(i);
        lTmpLong = of.getPeriod();

        if (lTmpLong == lPeriod)
        {
            //Ignore unknown variables
            if (of.getVariable() != geotop::input::UNKNOWN_VAR)
            {
                switch(of.getIntegrationType())
                {
                    case geotop::input::INS:
                        (lInstants->at(j))->append(of);
                        break;
                    case geotop::input::CUM:
                        (lCumulates->at(j))->append(of);
                        break;
                    case geotop::input::AVG:
                        (lAverages->at(j))->append(of);
                        break;
                    default:
                        break;
                }
            }

            i++;
        }
        else
        {
            lPeriod = lTmpLong;

            lInstants->push_back(new OutputFilesVector(lPeriod));
            lCumulates->push_back(new OutputFilesVector(lPeriod));
            lAverages->push_back(new OutputFilesVector(lPeriod));

            j++;
        }
    }

    //Pruning
    
    bool doneI = false;
    bool doneC = false;
    bool doneA = false;

    std::vector<OutputFilesVector*>::iterator itI = lInstants->end();
    std::vector<OutputFilesVector*>::iterator itC = lCumulates->end();
    std::vector<OutputFilesVector*>::iterator itA = lAverages->end();

    while (!doneI || !doneC || !doneA)
    {
        if (lInstants->empty())
            doneI = true;
        if (lCumulates->empty())
            doneC = true;
        if (lAverages->empty())
            doneA = true;

        if (!doneI)
        {
            if (itI != lInstants->begin())
                --itI;

            if ((*itI)->empty() == true)
            {
                itI = lInstants->erase(itI);
            }
            else if (itI == lInstants->begin())
                doneI = true;
        }

        if (!doneC)
        {
            if (itC != lCumulates->begin())
                --itC;

            if ((*itC)->empty() == true)
            {
                itC = lCumulates->erase(itC);
            }
            else if (itC == lCumulates->begin())
                doneC = true;
        }

        if (!doneA)
        {
            if (itA != lAverages->begin())
                --itA;

            if ((*itA)->empty() == true)
            {
                itA = lAverages->erase(itA);
            }
            else if (itA == lAverages->begin())
                doneA = true;
        }
    }

    if (lInstants->empty() == false)
        ofs.instants = lInstants;
    if (lCumulates->empty() == false)
        ofs.cumulates = lCumulates;
    if (lAverages->empty() == false)
        ofs.averages = lAverages;
    
}

void write_output_new(AllData* A)
{
    double lJDate = A->I->time; //seconds passed since the beginning of the simulation
    OutputFilesVector* ofv = NULL;
    size_t i, j;

    lJDate /= GTConst::secinday;
    lJDate += A->P->init_date;

    if (ofs.instants != NULL)
    {
        for (i = 0; i < ofs.instants->size(); i++)
        {
            ofv = ofs.instants->at(i);

            if (fabs(fmod(A->I->time, ofv->period())) < epsilon)
            {
                for (j = 0; j < ofv->size(); j++)
                {
                    std::cout << ofv->at(j).getFileName(convert_JDfrom0_dateeur12(lJDate)) << std::endl ;
                    geotop::input::OutputFile f = ofv->at(j);
                    printInstant(A, &f);
                }
            }
        }
    }

    if (ofs.cumulates != NULL)
    {
        for (i = 0; i < ofs.cumulates->size(); i++)
        {
            ofv = ofs.cumulates->at(i);

            if (fabs(fmod(A->I->time, ofv->period())) < epsilon)
            {
                for (j = 0; j < ofv->size(); j++)
                {
                    std::cout << ofv->at(j).getFileName(convert_JDfrom0_dateeur12(lJDate)) << std::endl ;
                }
            }
        }
    }

    if (ofs.averages != NULL)
    {
        for (i = 0; i < ofs.averages->size(); i++)
        {
            ofv = ofs.averages->at(i);

            if (fabs(fmod(A->I->time, ofv->period())) < epsilon)
            {
                for (j = 0; j < ofv->size(); j++)
                {
                    std::cout << ofv->at(j).getFileName(convert_JDfrom0_dateeur12(lJDate)) << std::endl ;
                }
            }
        }
    }

}

