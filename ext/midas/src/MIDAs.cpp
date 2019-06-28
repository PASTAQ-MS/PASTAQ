#include "MIDAs/MIDAs.h"
#include "ElementalComposition.hpp"

#if !defined WEB
#define WEB 1  // O==Not using WebInterface   1==Using WebInterface
#endif
#define cerr \
    if (!WEB) cerr
#define SWAP(a, b) \
    tempr = (a);   \
    (a) = (b);     \
    (b) = tempr
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
MIDAs::MIDAs(int charge_state, double fine_resolution, double coarse_resolution,
             double cgid_min_prob, int flag_method) {
    unsigned int i;
    struct Polynomial temp;

    if (flag_method != 1 && flag_method != 2 && flag_method != 3) {
        flag_method = 1;
    }
    if (cgid_min_prob > 1e-20) {
        cgid_min_prob = 1e-20;
    } else if (cgid_min_prob < 1e-300) {
        cgid_min_prob = 1e-300;
    }

    FLAG_METHOD = flag_method;  // flag of method used 1==P (Polynomial),2==FT
                                // (Foutier-Transform)
    FINE_RESOLUTION = fine_resolution;      // resolution of FGID
    COARSE_RESOLUTION = coarse_resolution;  // resolution of CGID
    MW_RESOLUTION = 1e-12;                  // resolution of molecular weight
    CHARGE_STATE = charge_state;            // charge of molecule
    MASS_ELECTRON = 0.00054857990946;       // electron mass
    MASS_HYDROGEN_ION = 1.0072764522;       // hydrogen ion mass
    MIN_PROB = cgid_min_prob;   // probability cut-off used for P_CGID
    FINE_MIN_PROB = 1e-16;      // probability cut-off used for P_FGID
    FINE_GRID = 1000;           // max number of mass kept during P_CGID
    FGID_TIME = CGID_TIME = 0;  // store algorithm time
    MONOISOTOPIC_MASS = 0;      // store monoisotopic mass

    if (flag_method == 1) {
        if (FINE_RESOLUTION > 0.1) {
            VECTOR_FGID_SIZE = unsigned(2000 / FINE_RESOLUTION + 0.5);
        } else {
            VECTOR_FGID_SIZE = unsigned(500 / FINE_RESOLUTION + 0.5);
        }

        for (i = 0; i <= VECTOR_FGID_SIZE; i++) {
            temp.prob = 0;
            temp.power = 0;
            FGID_Polynomial.push_back(
                temp);  // used for Fine-grained polynomial multiplication
        }

        VECTOR_CGID_SIZE = unsigned(2000);
        for (i = 0; i <= VECTOR_CGID_SIZE; i++) {
            temp.prob = 0;
            temp.power = 0;
            CGID_Polynomial.push_back(
                temp);  // used for Coarse-grained polynomial multiplication
        }
    } else {
        VECTOR_FGID_SIZE = VECTOR_CGID_SIZE = 0;
    }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*Initializing Elemental Composition*/
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void MIDAs::Initialize_Elemental_Composition(
    string substance, string mc, string nt, string ct, int flag_substance,
    vector<struct Composition> user_element_composition) {
    int i;

    User_Element_Composition = user_element_composition;

    // Example Element Composition Class
    // ElementComposition example( "Cysteine PTM","NT","CT" ); //Cesteine is not
    // modified, N-terminal, C-terminal
    ElementComposition ec(&mc[0], &nt[0], &ct[0], flag_substance);

    if (User_Element_Composition.size() > 0) {
        // Updating user elemental composition
        cerr << "ID:Start user specify elemental composition" << endl;
        for (i = 0; i < int(User_Element_Composition.size()); i++) {
            ec.UpdateElement(User_Element_Composition[i]);
        }
        cerr << "ID:Finished user specify elemental composition" << endl;
    }

    cerr << "ID:Start initializing elemental composition. Substance="
         << substance << endl;
    ec.CalcComposition(&substance[0], Element_Composition);
    cerr << "ID:Finished initializing elemental composition" << endl;

    cerr << "ID:Start correcting charge state" << endl;
    add_charge_state(Element_Composition, CHARGE_STATE);
    cerr << "ID:Finished correcting charge state" << endl;

    cerr << "ID:Start sorting element composition" << endl;
    sort_element_composition(Element_Composition);
    cerr << "ID:Finished sorting element composition" << endl;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*Setting resolution*/
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void MIDAs::set_resolution(double mw_mono) {
    double coarse_resolution, fine_resolution;
    cerr.precision(12);
    cerr << "ID:---------------------------------------------------------------"
            "-------------------------"
         << endl;
    cerr << "ID:Lightest mass = " << mw_mono << endl;
    cerr << "ID:CGID resolution = " << COARSE_RESOLUTION << endl;
    cerr << "ID:CGID minimum probability = " << MIN_PROB << endl;
    cerr << "ID:FGID resolution = " << FINE_RESOLUTION << endl;
    cerr << "ID:FGID minimum probability = " << FINE_MIN_PROB << endl;
    cerr << "ID:MW resolution = " << MW_RESOLUTION << endl;

    fine_resolution = FINE_RESOLUTION;
    if (COARSE_RESOLUTION > 10) {
        COARSE_RESOLUTION = 10;
    }
    coarse_resolution = COARSE_RESOLUTION;

    // Resolution used to merge coarse-ID
    if (coarse_resolution <= 1.0) {
        MERGE_COARSE_RESOLUTION = (1.0 - 0.022);
    } else {
        MERGE_COARSE_RESOLUTION = (coarse_resolution)-0.1;
    }

    // Resolution to compute coarse-ID
    if (coarse_resolution <= 1.0 || coarse_resolution >= 1.0) {
        coarse_resolution = 0.9;
    }

    if (fine_resolution >= 1.0) {
        MERGE_FINE_RESOLUTION = (1.0) - 0.022;
        fine_resolution = 0.9;
    } else if (fine_resolution <= 1e-4 && mw_mono < 1e5) {
        fine_resolution = 1e-4;
        MERGE_FINE_RESOLUTION = fine_resolution;
    } else if (fine_resolution <= 1e-3 && mw_mono < 1e6) {
        fine_resolution = 1e-3;
        MERGE_FINE_RESOLUTION = fine_resolution;
    } else if (fine_resolution <= 1e-2 && mw_mono < 2e6) {
        fine_resolution = 1e-2;
        MERGE_FINE_RESOLUTION = fine_resolution;
    } else {
        MERGE_FINE_RESOLUTION = fine_resolution;
        fine_resolution = 1e-2;
    }

    FINE_RESOLUTION = fine_resolution / 2.0;
    COARSE_RESOLUTION = coarse_resolution / 2.0;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*Returns computed fine-grained isotopic distribution*/
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
vector<struct Isotopic_Distribution>
MIDAs::Fine_Grained_Isotopic_Distribution() {
    int i;
    double np, mw_mono;
    vector<struct Polynomial> F_Polynomial, T_fgid;
    struct Isotopic_Distribution temp;
    struct timeval start_n, end_n;

    cerr.precision(16);
    cout.precision(16);
    cerr << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
            "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
         << endl;
    MONOISOTOPIC_MASS = mw_mono =
        monoisotopic_molecular_mass(Element_Composition);
    cerr << "ID:Finished computing lightest mass" << endl;
    set_resolution(mw_mono);
    cerr << "ID:Finished setting resolution" << endl;
    cerr << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
            "%%%%%%%%%%%%%%%%%%%%%%%%%%"
         << endl;

    /*Computing Coarse-Grain ID.*/
    if (FLAG_METHOD == 1)  // Multiplying Polynomial
    {
        /*Computing Fine-Grained Isotopic Distribution*/
        cerr << "ID:Started polynomial multiplication method-1." << endl;
        cerr << "ID:Computing****Fine-Grained Isotopic Distribution (FGID).****"
             << endl;
        gettimeofday(&start_n, NULL);
        T_fgid = multiply_fine_polynomial(Element_Composition);
        T_fgid = Merge_Fine_Polynomial(T_fgid);
        cerr << "ID:Finished computing ****FGID****" << endl;
        gettimeofday(&end_n, NULL);
        FGID_TIME = (end_n.tv_sec + end_n.tv_usec / 1000000.0) -
                    (start_n.tv_sec + start_n.tv_usec / 1000000.0);
        cerr << "ID:***********************************************************"
                "*****************************"
             << endl;
        cerr << "ID:Elapsed time (seconds) FGID = " << FGID_TIME << endl;
        cerr << "ID:***********************************************************"
                "*****************************"
             << endl;

        /*Computing sum of probability from computed ID np*/
        np = 0;
        for (i = 0; i < int(T_fgid.size()); i++) {
            if (CHARGE_STATE != 0) {
                T_fgid[i].power =
                    (T_fgid[i].power * MW_RESOLUTION) / fabs(CHARGE_STATE);
            } else {
                T_fgid[i].power = (T_fgid[i].power * MW_RESOLUTION);
            }
            np = np + T_fgid[i].prob;

            temp.mw = T_fgid[i].power;
            temp.prob = T_fgid[i].prob;
            FINE_GRAINED_ISOTOPIC_DISTRIBUTION.push_back(temp);
        }

        // normalize probability
        for (i = 0; i < int(FINE_GRAINED_ISOTOPIC_DISTRIBUTION.size()); i++) {
            FINE_GRAINED_ISOTOPIC_DISTRIBUTION[i].prob =
                FINE_GRAINED_ISOTOPIC_DISTRIBUTION[i].prob / np;
        }

        cerr << "ID:Sum of probabilities FGID = " << np << endl;
    } else if (FLAG_METHOD == 2)  // Fourier Transform Method-1
    {
        /********Fine-Grained IS****************/
        cerr << "ID:Started Fourier Transform-1." << endl;
        cerr
            << "ID:Computing****Fine-Grained Isotopic Distribution (FGID).*****"
            << endl;
        cerr << "ID:-----------------------------------------------------------"
                "-----------------------------"
             << endl;
        gettimeofday(&start_n, NULL);
        FT_Fine_Grained_ID(Element_Composition, T_fgid, FINE_RESOLUTION);
        cerr << "ID:Finished computing****FGID****" << endl;
        gettimeofday(&end_n, NULL);
        FGID_TIME = (end_n.tv_sec + end_n.tv_usec / 1000000.0) -
                    (start_n.tv_sec + start_n.tv_usec / 1000000.0);
        cerr << "ID:***********************************************************"
                "*****************************"
             << endl;
        cerr << "ID:Elapsed time (seconds) FGID = " << FGID_TIME << endl;
        cerr << "ID:***********************************************************"
                "*****************************"
             << endl;

        for (i = 0; i < int(T_fgid.size()) - 1; i++) {
            if (CHARGE_STATE != 0) {
                T_fgid[i].power = T_fgid[i].power / fabs(CHARGE_STATE);
            }
            temp.mw = T_fgid[i].power;
            temp.prob = T_fgid[i].prob;
            FINE_GRAINED_ISOTOPIC_DISTRIBUTION.push_back(temp);
        }
    }
    return FINE_GRAINED_ISOTOPIC_DISTRIBUTION;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*Returns fine-grained isotopic distribution*/
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
vector<struct Isotopic_Distribution>
MIDAs::Coarse_Grained_Isotopic_Distribution() {
    int i;
    double np, sum_mw, sum_p, mw_mono;
    vector<struct Polynomial> F_Polynomial;
    struct Isotopic_Distribution temp;
    struct timeval start_n, end_n;

    cerr.precision(16);
    cout.precision(16);
    /*Reading Element Composition*/
    cerr << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
            "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
         << endl;
    MONOISOTOPIC_MASS = mw_mono =
        monoisotopic_molecular_mass(Element_Composition);
    cerr << "ID:Finished computing lightest mass" << endl;
    set_resolution(mw_mono);
    cerr << "ID:Finished setting resolution" << endl;
    cerr << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
            "%%%%%%%%%%%%%%%%%%%%%%%%%%"
         << endl;

    /*Computing Coarse-Grain ID.*/
    if (FLAG_METHOD == 1)  // Multiplying Polynomial
    {
        cerr << "ID:Started polynomial multiplication method-1." << endl;
        cerr << "ID:Computing****Coarse-Grained Isotopic Distribution "
                "(CGID).****"
             << endl;
        gettimeofday(&start_n, NULL);
        F_Polynomial = multiply_polynomial(Element_Composition);
        F_Polynomial = Merge_Coarse_Polynomial(F_Polynomial);
        cerr << "ID:Finished computing ****CGID****" << endl;
        gettimeofday(&end_n, NULL);
        CGID_TIME = (end_n.tv_sec + end_n.tv_usec / 1000000.0) -
                    (start_n.tv_sec + start_n.tv_usec / 1000000.0);
        cerr << "ID:***********************************************************"
                "*****************************"
             << endl;
        cerr << "ID:Elapsed time (seconds) CGID = " << CGID_TIME << endl;
        cerr << "ID:***********************************************************"
                "*****************************"
             << endl;

        /*Computing normalizing probability factor np*/
        np = 0;
        for (i = 0; i < int(F_Polynomial.size()); i++) {
            np = np + F_Polynomial[i].prob;
        }
        cerr << "ID:Sum of probabilities CGID = " << np << endl;

        sum_mw = sum_p = 0.0;
        for (i = 0; i < int(F_Polynomial.size()); i++) {
            temp.mw = F_Polynomial[i].power * MW_RESOLUTION;
            temp.prob = F_Polynomial[i].prob / np;

            sum_mw = sum_mw + temp.mw * temp.prob;
            sum_p = sum_p + temp.prob;

            COARSE_GRAINED_ISOTOPIC_DISTRIBUTION.push_back(temp);
        }
        cerr << "ID:Sum of probabilities CGID = " << sum_p << endl;
    } else if (FLAG_METHOD == 2)  // Fourier Transform Method-1
    {
        cerr << "ID:Started Fourier Transform-1." << endl;
        cerr << "ID:Computing****Coarse-Grained Isotopic Distribution "
                "(CGID).****"
             << endl;
        gettimeofday(&start_n, NULL);
        COARSE_GRAINED_ISOTOPIC_DISTRIBUTION.clear();
        FT_Coarse_Grained_ID(Element_Composition,
                             COARSE_GRAINED_ISOTOPIC_DISTRIBUTION,
                             COARSE_RESOLUTION);
        cerr << "ID:Finished computing ****CGID****" << endl;
        gettimeofday(&end_n, NULL);
        CGID_TIME = (end_n.tv_sec + end_n.tv_usec / 1000000.0) -
                    (start_n.tv_sec + start_n.tv_usec / 1000000.0);
        cerr << "ID:***********************************************************"
                "*****************************"
             << endl;
        cerr << "ID:Elapsed time (seconds) CGID = " << CGID_TIME << endl;
        cerr << "ID:***********************************************************"
                "*****************************"
             << endl;
        cerr << "ID:-----------------------------------------------------------"
                "-----------------------------"
             << endl;
    } else if (FLAG_METHOD == 3)  // JFC method FTT-Method-2
    {
        cerr << "ID:Started Fourier Transform Method-2." << endl;
        cerr << "ID:Computing****Coarse-Grained Isotopic Distribution "
                "(CGID).****"
             << endl;
        gettimeofday(&start_n, NULL);
        COARSE_GRAINED_ISOTOPIC_DISTRIBUTION.clear();
        FT_JFC_ID(Element_Composition);
        cerr << "ID:Finished computing ****CGID****" << endl;
        gettimeofday(&end_n, NULL);
        CGID_TIME = (end_n.tv_sec + end_n.tv_usec / 1000000.0) -
                    (start_n.tv_sec + start_n.tv_usec / 1000000.0);
        cerr << "ID:***********************************************************"
                "*****************************"
             << endl;
        cerr << "ID:Elapsed time (seconds) CGID = " << CGID_TIME << endl;
        cerr << "ID:***********************************************************"
                "*****************************"
             << endl;
        cerr << "ID:-----------------------------------------------------------"
                "-----------------------------"
             << endl;
    }
    cerr << "ID:---------------------------------------------------------------"
            "-------------------------"
         << endl;

    return COARSE_GRAINED_ISOTOPIC_DISTRIBUTION;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*Printing fine structure isotopic distirbution*/
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void MIDAs::Print_Fine_Structure_Isotopic_Distribution(char *filename) {
    int i, count;
    double mw_mono, sum_p = 0.0;
    char outfile[256];
    FILE *out = NULL;

    mw_mono = monoisotopic_molecular_mass(Element_Composition);

#if !WEB
    strcpy(outfile, filename);
    strcat(outfile, "_Fine_Grained_Isotopic_Distribution");
    out = fopen(outfile, "w+");

    count = 0;
    for (i = 0; i < int(FINE_GRAINED_ISOTOPIC_DISTRIBUTION.size()); i++) {
        if (FINE_GRAINED_ISOTOPIC_DISTRIBUTION[i].mw != 0 &&
            FINE_GRAINED_ISOTOPIC_DISTRIBUTION[i].mw >= (mw_mono - 0.01)) {
            count = count + 1;
            sum_p = sum_p + FINE_GRAINED_ISOTOPIC_DISTRIBUTION[i].prob;

            fprintf(out, "%.16f\t%.16e\n",
                    FINE_GRAINED_ISOTOPIC_DISTRIBUTION[i].mw,
                    FINE_GRAINED_ISOTOPIC_DISTRIBUTION[i].prob);
        }
    }
    fclose(out);
    cerr << "ID:---------------------------------------------------------------"
            "-------------------------"
         << endl;
    cerr << "ID:Number of terms Fine-grained ID = " << count << endl;
    cerr << "ID:---------------------------------------------------------------"
            "-------------------------"
         << endl;
    // cerr<<"ID:Sum of probabilities for fine structure isotopic distribution =
    // "<<sum_p<<endl;
#else
    out = stdout;
    strcpy(outfile, "Web output");

    fprintf(out, "\\Lightest_theoretical_molecular_mass\n");
    fprintf(out, "%.16lf\n", monoisotopic_molecular_mass(Element_Composition));

    fprintf(out, "\\Average_theoretical_mass\n");
    fprintf(out, "%.16lf\n", average_molecular_mass());

    if (COARSE_GRAINED_ISOTOPIC_DISTRIBUTION.size() > 0) {
        fprintf(out, "\\Average_computed_molecular_mass\n");
        fprintf(out, "%.16lf\n",
                average_molecular_mass_distribution(
                    COARSE_GRAINED_ISOTOPIC_DISTRIBUTION));
    } else if (FINE_GRAINED_ISOTOPIC_DISTRIBUTION.size() > 0) {
        fprintf(out, "\\Average_computed_molecular_mass\n");
        fprintf(out, "%.16lf\n",
                average_molecular_mass_distribution(
                    FINE_GRAINED_ISOTOPIC_DISTRIBUTION));
    }

    fprintf(out, "\\Theoretical_standard_deviation\n");
    fprintf(out, "%.16lf\n", sqrt(variance_molecular_mass()));

    if (COARSE_GRAINED_ISOTOPIC_DISTRIBUTION.size() > 0) {
        fprintf(out, "\\Computed_standard_deviation\n");
        fprintf(out, "%.16lf\n",
                sqrt(variance_molecular_mass_distribution(
                    COARSE_GRAINED_ISOTOPIC_DISTRIBUTION)));
    } else if (FINE_GRAINED_ISOTOPIC_DISTRIBUTION.size() > 0) {
        fprintf(out, "\\Computed_standard_deviation\n");
        fprintf(out, "%.16lf\n",
                sqrt(variance_molecular_mass_distribution(
                    FINE_GRAINED_ISOTOPIC_DISTRIBUTION)));
    }

    fprintf(
        out,
        "\\Fine_Grained_Isotopic_Distribution:(Molecular_Mass,Probability)\n");

    for (i = 0; i < int(FINE_GRAINED_ISOTOPIC_DISTRIBUTION.size()); i++) {
        if (FINE_GRAINED_ISOTOPIC_DISTRIBUTION[i].mw != 0 &&
            FINE_GRAINED_ISOTOPIC_DISTRIBUTION[i].mw >= (mw_mono - 0.01)) {
            if (CHARGE_STATE != 0) {
                fprintf(out, "%.16f\t%.16e\n",
                        FINE_GRAINED_ISOTOPIC_DISTRIBUTION[i].mw,
                        FINE_GRAINED_ISOTOPIC_DISTRIBUTION[i].prob);
            } else {
                fprintf(out, "%.16f\t%.16e\n",
                        FINE_GRAINED_ISOTOPIC_DISTRIBUTION[i].mw,
                        FINE_GRAINED_ISOTOPIC_DISTRIBUTION[i].prob);
            }
        }
    }

    fprintf(out,
            "\\Coarse_Grained_Isotopic_Distribution:(Molecular_Mass,"
            "Probability)\n");

    for (i = 0; i < int(COARSE_GRAINED_ISOTOPIC_DISTRIBUTION.size()); i++) {
        if (COARSE_GRAINED_ISOTOPIC_DISTRIBUTION[i].mw >= (mw_mono - 0.01)) {
            if (CHARGE_STATE != 0) {
                fprintf(out, "%.16f  %.16e\n",
                        COARSE_GRAINED_ISOTOPIC_DISTRIBUTION[i].mw /
                            fabs(CHARGE_STATE),
                        COARSE_GRAINED_ISOTOPIC_DISTRIBUTION[i].prob);
            } else {
                fprintf(out, "%.16f  %.16e\n",
                        COARSE_GRAINED_ISOTOPIC_DISTRIBUTION[i].mw,
                        COARSE_GRAINED_ISOTOPIC_DISTRIBUTION[i].prob);
            }
        }
    }
    fclose(out);

#endif
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*Printing average structure isotopic distribution*/
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void MIDAs::Print_Coarse_Structure_Isotopic_Distribution(char *filename) {
    int i, count;
    double mw_mono;
    char outfile[256];
    FILE *out = NULL;

    mw_mono = monoisotopic_molecular_mass(Element_Composition);

#if !WEB
    strcpy(outfile, filename);
    strcat(outfile, "_Coarse_Grained_Isotopic_Distribution");
    out = fopen(outfile, "w+");
    count = 0;
    for (i = 0; i < int(COARSE_GRAINED_ISOTOPIC_DISTRIBUTION.size()); i++) {
        if (COARSE_GRAINED_ISOTOPIC_DISTRIBUTION[i].mw >= (mw_mono - 0.01)) {
            count = count + 1;
            if (CHARGE_STATE != 0) {
                fprintf(out, "%.16f  %.16e\n",
                        COARSE_GRAINED_ISOTOPIC_DISTRIBUTION[i].mw /
                            fabs(CHARGE_STATE),
                        COARSE_GRAINED_ISOTOPIC_DISTRIBUTION[i].prob);
            } else {
                fprintf(out, "%.16f  %.16e\n",
                        COARSE_GRAINED_ISOTOPIC_DISTRIBUTION[i].mw,
                        COARSE_GRAINED_ISOTOPIC_DISTRIBUTION[i].prob);
            }
        }
    }
    fclose(out);
    cerr << "ID:---------------------------------------------------------------"
            "-------------------------"
         << endl;
    cerr << "ID:Number of terms Coarse-grained ID = " << count << endl;
    cerr << "ID:---------------------------------------------------------------"
            "-------------------------"
         << endl;
#endif
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*Computing the fine isotopic distribution*/
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
vector<struct Polynomial> MIDAs::multiply_fine_polynomial(
    vector<vector<struct Composition> > &Element_Composition) {
    int i, j, k, n, max_isotope;
    int nc = 10.0;  // 7 gp = 1e-11, 8 gp=1e-14, 9 gp 1e-16
    int nc_add = 1;
    int nc_add_value = 1.0;
    int n_atoms = 200;
    double max_polynomial_size =
        log(1e13);  // max munber of term allowed to used polynomial method
    double n_polynomial_terms = 0;
    struct Polynomial temp_Polynomial;
    vector<struct Polynomial> T_Polynomial;
    vector<vector<struct Composition> > element_composition;
    double min_prob;
    double prob;

    cout.precision(16);

    // Probability cut-off
    min_prob = FINE_MIN_PROB;

    // Count number of unique elements and transform probabilities to log_base
    max_isotope = n = 0;
    for (k = 0; k < int(Element_Composition.size()); k++) {
        if (Element_Composition[k].size() > 0) {
            n = n + 1;
        }
        if (Element_Composition[k].size() > 10) {
            max_isotope = 1;
        }
    }
    vector<vector<struct Polynomial> > F_Polynomial(n);

    cerr << "ID:Polynomial expansion in progress" << endl;
    cerr << "ID:---------------------------------------------------------------"
            "-------------------------"
         << endl;

    if (max_isotope == 0) {
        for (k = 0; k < n; k++) {
            T_Polynomial.clear();
            // cout<<"ID::Fine Structure "<<k<<endl;
            if (Element_Composition[k].size() ==
                1)  // Handles the case of two isotopic elements
            {
                int n1;
                // computing mean and standard deviation
                n1 = int(Element_Composition[k][0].atoms *
                         Element_Composition[k][0].prob);
                // ns1=int(ceil(nc*sqrt(Element_Composition[k][0].atoms*Element_Composition[k][0].prob*(1.0-Element_Composition[k][0].prob))));
                // cout<<n1<<"\t"<<ns1<<"\t"<<Element_Composition[k][0].prob<<endl;

                prob = FACTR_LN(Element_Composition[k][0].atoms) -
                       FACTR_LN(n1) + n1 * (Element_Composition[k][0].log_prob);
                prob = exp(prob);

                temp_Polynomial.power = n1 * Element_Composition[k][0].power;
                temp_Polynomial.prob = prob;

                F_Polynomial[k].push_back(temp_Polynomial);

                cerr << "ID:Number of isotopes 1. Polynomial size = "
                     << F_Polynomial[k].size() << endl;
            } else if (Element_Composition[k].size() ==
                       2)  // Handles the case of two isotopic elements
            {
                if (Element_Composition[k][0].atoms < n_atoms) {
                    nc_add = 10;
                } else {
                    nc_add = nc_add_value;
                }

                int n1, n2, ns1, ns2;
                // computing mean and standard deviation
                n1 = int(Element_Composition[k][0].atoms *
                         Element_Composition[k][0].prob);
                n2 = int(Element_Composition[k][0].atoms *
                         Element_Composition[k][1].prob);

                ns1 = int(
                    ceil(nc_add +
                         nc * sqrt(Element_Composition[k][0].atoms *
                                   Element_Composition[k][0].prob *
                                   (1.0 - Element_Composition[k][0].prob))));
                ns2 = int(
                    ceil(nc_add +
                         nc * sqrt(Element_Composition[k][1].atoms *
                                   Element_Composition[k][1].prob *
                                   (1.0 - Element_Composition[k][1].prob))));
                n2 = n2 + ns2;

                n_polynomial_terms = log(n1 + ns1) + log(n2) + log(pow(2, 2));
                if (n_polynomial_terms >
                    max_polynomial_size)  // checking polynomial size
                {
                    element_composition.push_back(Element_Composition[k]);
                    FT_Fine_Grained_ID(element_composition, T_Polynomial,
                                       FINE_RESOLUTION);
                    for (i = 0; i < int(T_Polynomial.size()); i++) {
                        if (T_Polynomial[i].power > 0) {
                            T_Polynomial[i].power =
                                T_Polynomial[i].power / MW_RESOLUTION;
                            F_Polynomial[k].push_back(T_Polynomial[i]);
                        }
                    }
                    element_composition.clear();
                    T_Polynomial.clear();
                } else {
                    for (i = (n1 + ns1); i >= (n1 - ns1) && i >= 0; i--) {
                        j = Element_Composition[k][0].atoms - i;

                        if (j >= 0) {
                            prob = FACTR_LN(Element_Composition[k][0].atoms) -
                                   FACTR_LN(i) - FACTR_LN(j) +
                                   i * (Element_Composition[k][0].log_prob) +
                                   j * (Element_Composition[k][1].log_prob);

                            prob = exp(prob);

                            if (prob >= min_prob) {
                                temp_Polynomial.power =
                                    i * Element_Composition[k][0].power +
                                    j * Element_Composition[k][1].power;
                                temp_Polynomial.prob = prob;
                                F_Polynomial[k].push_back(temp_Polynomial);
                            }
                        }
                    }
                }
                cerr << "ID:" << n1 + ns1 << "," << n2 << endl;
                cerr << "ID:Number of isotopes 2. Polynomial size = "
                     << F_Polynomial[k].size() << endl;
            } else if (Element_Composition[k].size() ==
                       3)  // Handles the case of three isotopic elements
            {
                if (Element_Composition[k][0].atoms < n_atoms) {
                    nc_add = 10;
                } else {
                    nc_add = nc_add_value;
                }

                int n1, n2, n3, ns1, ns2, ns3, l;

                // computing mean and standard deviation
                n1 = int(Element_Composition[k][0].atoms *
                         Element_Composition[k][0].prob);
                n2 = int(Element_Composition[k][0].atoms *
                         Element_Composition[k][1].prob);
                n3 = int(Element_Composition[k][0].atoms *
                         Element_Composition[k][2].prob);

                ns1 = int(
                    ceil(nc_add +
                         nc * sqrt(Element_Composition[k][0].atoms *
                                   Element_Composition[k][0].prob *
                                   (1.0 - Element_Composition[k][0].prob))));
                ns2 = int(
                    ceil(nc_add +
                         nc * sqrt(Element_Composition[k][1].atoms *
                                   Element_Composition[k][1].prob *
                                   (1.0 - Element_Composition[k][1].prob))));
                ns3 = int(
                    ceil(nc_add +
                         nc * sqrt(Element_Composition[k][2].atoms *
                                   Element_Composition[k][2].prob *
                                   (1.0 - Element_Composition[k][2].prob))));
                n3 = n3 + ns3;

                n_polynomial_terms =
                    log(n1 + ns1) + log(n2 + ns2) + log(n3) + log(pow(2, 3));
                if (n_polynomial_terms >
                    max_polynomial_size)  // checking polynomial size
                {
                    element_composition.push_back(Element_Composition[k]);
                    FT_Fine_Grained_ID(element_composition, T_Polynomial,
                                       FINE_RESOLUTION);
                    for (i = 0; i < int(T_Polynomial.size()); i++) {
                        if (T_Polynomial[i].power > 0) {
                            T_Polynomial[i].power =
                                T_Polynomial[i].power / MW_RESOLUTION;
                            F_Polynomial[k].push_back(T_Polynomial[i]);
                        }
                    }
                    element_composition.clear();
                    T_Polynomial.clear();
                } else {
                    for (i = n1 + ns1; i >= n1 - ns1 && i >= 0; i--) {
                        for (j = n2 + ns2; j >= n2 - ns2 && j >= 0; j--) {
                            l = Element_Composition[k][0].atoms - (i + j);

                            if (l >= 0 && (l <= n3)) {
                                prob =
                                    FACTR_LN(Element_Composition[k][0].atoms) -
                                    FACTR_LN(i) - FACTR_LN(j) - FACTR_LN(l) +
                                    i * (Element_Composition[k][0].log_prob) +
                                    j * (Element_Composition[k][1].log_prob) +
                                    l * (Element_Composition[k][2].log_prob);

                                prob = exp(prob);

                                if (prob >= min_prob) {
                                    temp_Polynomial.power =
                                        i * Element_Composition[k][0].power +
                                        j * Element_Composition[k][1].power +
                                        l * Element_Composition[k][2].power;
                                    temp_Polynomial.prob = prob;

                                    F_Polynomial[k].push_back(temp_Polynomial);
                                }
                            }
                        }
                    }
                }
                cerr << "ID:" << n1 + ns1 << "," << n2 + ns2 << "," << n3
                     << endl;
                cerr << "ID:Number of isotopes 3. Polynomial size = "
                     << F_Polynomial[k].size() << endl;
            } else if (Element_Composition[k].size() ==
                       4)  // Handles the case of four isotopic elements
            {
                if (Element_Composition[k][0].atoms < n_atoms) {
                    nc_add = 10;
                } else {
                    nc_add = nc_add_value;
                }

                int n1, n2, n3, n4, ns1, ns2, ns3, ns4, l, m;

                // computing mean and standard deviation
                n1 = int(Element_Composition[k][0].atoms *
                         Element_Composition[k][0].prob);
                n2 = int(Element_Composition[k][0].atoms *
                         Element_Composition[k][1].prob);
                n3 = int(Element_Composition[k][0].atoms *
                         Element_Composition[k][2].prob);
                n4 = int(Element_Composition[k][0].atoms *
                         Element_Composition[k][3].prob);

                ns1 = int(
                    ceil(nc_add +
                         nc * sqrt(Element_Composition[k][0].atoms *
                                   Element_Composition[k][0].prob *
                                   (1.0 - Element_Composition[k][0].prob))));
                ns2 = int(
                    ceil(nc_add +
                         nc * sqrt(Element_Composition[k][1].atoms *
                                   Element_Composition[k][1].prob *
                                   (1.0 - Element_Composition[k][1].prob))));
                ns3 = int(
                    ceil(nc_add +
                         nc * sqrt(Element_Composition[k][2].atoms *
                                   Element_Composition[k][2].prob *
                                   (1.0 - Element_Composition[k][2].prob))));
                ns4 = int(
                    ceil(nc_add +
                         nc * sqrt(Element_Composition[k][3].atoms *
                                   Element_Composition[k][3].prob *
                                   (1.0 - Element_Composition[k][3].prob))));

                n4 = n4 + ns4;
                n_polynomial_terms = log(n1 + ns1) + log(n2 + ns2) +
                                     log(n3 + ns3) + log(n4) + log(pow(2, 4));
                if (n_polynomial_terms >
                    max_polynomial_size)  // checking polynomial size
                {
                    element_composition.push_back(Element_Composition[k]);
                    FT_Fine_Grained_ID(element_composition, T_Polynomial,
                                       FINE_RESOLUTION);
                    for (i = 0; i < int(T_Polynomial.size()); i++) {
                        if (T_Polynomial[i].power > 0) {
                            T_Polynomial[i].power =
                                T_Polynomial[i].power / MW_RESOLUTION;
                            F_Polynomial[k].push_back(T_Polynomial[i]);
                        }
                    }
                    element_composition.clear();
                    T_Polynomial.clear();
                } else {
                    for (i = n1 + ns1; i >= n1 - ns1 && i >= 0; i--) {
                        for (j = n2 + ns2; j >= n2 - ns2 && j >= 0; j--) {
                            for (m = n3 + ns3; m >= n3 - ns3 && m >= 0; m--) {
                                l = Element_Composition[k][0].atoms -
                                    (i + j + m);
                                if ((l >= 0) && (l <= n4)) {
                                    prob =
                                        FACTR_LN(
                                            Element_Composition[k][0].atoms) -
                                        FACTR_LN(i) - FACTR_LN(j) -
                                        FACTR_LN(m) - FACTR_LN(l) +
                                        i * (Element_Composition[k][0]
                                                 .log_prob) +
                                        j * (Element_Composition[k][1]
                                                 .log_prob) +
                                        m * (Element_Composition[k][2]
                                                 .log_prob) +
                                        l * (Element_Composition[k][3]
                                                 .log_prob);
                                    prob = exp(prob);

                                    if (prob >= min_prob) {
                                        temp_Polynomial.power =
                                            i * Element_Composition[k][0]
                                                    .power +
                                            j * Element_Composition[k][1]
                                                    .power +
                                            m * Element_Composition[k][2]
                                                    .power +
                                            l * Element_Composition[k][3].power;
                                        temp_Polynomial.prob = prob;

                                        F_Polynomial[k].push_back(
                                            temp_Polynomial);
                                    }
                                }
                            }
                        }
                    }
                }
                cerr << "ID:" << n1 + ns1 << "," << n2 + ns2 << "," << n3 + ns3
                     << "," << n4 << endl;
                cerr << "ID:Number of isotopes 4. Polynomial size = "
                     << F_Polynomial[k].size() << endl;
            } else if (Element_Composition[k].size() ==
                       5)  // Handles the case of five isotopic elements
            {
                if (Element_Composition[k][0].atoms < n_atoms) {
                    nc_add = 10;
                } else {
                    nc_add = nc_add_value;
                }

                int n1, n2, n3, n4, n5, ns1, ns2, ns3, ns4, ns5, l, m, n;

                // computing mean and standard deviation
                n1 = int(Element_Composition[k][0].atoms *
                         Element_Composition[k][0].prob);
                n2 = int(Element_Composition[k][0].atoms *
                         Element_Composition[k][1].prob);
                n3 = int(Element_Composition[k][0].atoms *
                         Element_Composition[k][2].prob);
                n4 = int(Element_Composition[k][0].atoms *
                         Element_Composition[k][3].prob);
                n5 = int(Element_Composition[k][0].atoms *
                         Element_Composition[k][4].prob);

                ns1 = int(
                    ceil(nc_add +
                         nc * sqrt(Element_Composition[k][0].atoms *
                                   Element_Composition[k][0].prob *
                                   (1.0 - Element_Composition[k][0].prob))));
                ns2 = int(
                    ceil(nc_add +
                         nc * sqrt(Element_Composition[k][1].atoms *
                                   Element_Composition[k][1].prob *
                                   (1.0 - Element_Composition[k][1].prob))));
                ns3 = int(
                    ceil(nc_add +
                         nc * sqrt(Element_Composition[k][2].atoms *
                                   Element_Composition[k][2].prob *
                                   (1.0 - Element_Composition[k][2].prob))));
                ns4 = int(
                    ceil(nc_add +
                         nc * sqrt(Element_Composition[k][3].atoms *
                                   Element_Composition[k][3].prob *
                                   (1.0 - Element_Composition[k][3].prob))));
                ns5 = int(
                    ceil(nc_add +
                         nc * sqrt(Element_Composition[k][4].atoms *
                                   Element_Composition[k][4].prob *
                                   (1.0 - Element_Composition[k][4].prob))));

                n5 = n5 + ns5;
                n_polynomial_terms = log(n1 + ns1) + log(n2 + ns2) +
                                     log(n3 + ns3) + log(n4 + ns4) + log(n5) +
                                     log(pow(2, 5));
                if (n_polynomial_terms >
                    max_polynomial_size)  // checking polynomial size
                {
                    element_composition.push_back(Element_Composition[k]);
                    FT_Fine_Grained_ID(element_composition, T_Polynomial,
                                       FINE_RESOLUTION);
                    for (i = 0; i < int(T_Polynomial.size()); i++) {
                        if (T_Polynomial[i].power > 0) {
                            T_Polynomial[i].power =
                                T_Polynomial[i].power / MW_RESOLUTION;
                            F_Polynomial[k].push_back(T_Polynomial[i]);
                        }
                    }
                    element_composition.clear();
                    T_Polynomial.clear();
                } else {
                    for (i = n1 + ns1; i >= n1 - ns1 && i >= 0; i--) {
                        for (j = n2 + ns2; j >= n2 - ns2 && j >= 0; j--) {
                            for (m = n3 + ns3; m >= n3 - ns3 && m >= 0; m--) {
                                for (n = n4 + ns4; n >= n4 - ns4 && n >= 0;
                                     n--) {
                                    l = Element_Composition[k][0].atoms -
                                        (i + j + m + n);
                                    if ((l >= 0) && (l <= n5)) {
                                        prob =
                                            FACTR_LN(Element_Composition[k][0]
                                                         .atoms) -
                                            FACTR_LN(i) - FACTR_LN(j) -
                                            FACTR_LN(m) - FACTR_LN(n) -
                                            FACTR_LN(l) +
                                            i * (Element_Composition[k][0]
                                                     .log_prob) +
                                            j * (Element_Composition[k][1]
                                                     .log_prob) +
                                            m * (Element_Composition[k][2]
                                                     .log_prob) +
                                            n * (Element_Composition[k][3]
                                                     .log_prob) +
                                            l * (Element_Composition[k][4]
                                                     .log_prob);

                                        prob = exp(prob);

                                        if (prob >= min_prob) {
                                            temp_Polynomial.power =
                                                i * Element_Composition[k][0]
                                                        .power +
                                                j * Element_Composition[k][1]
                                                        .power +
                                                m * Element_Composition[k][2]
                                                        .power +
                                                n * Element_Composition[k][3]
                                                        .power +
                                                l * Element_Composition[k][4]
                                                        .power;
                                            temp_Polynomial.prob = prob;

                                            F_Polynomial[k].push_back(
                                                temp_Polynomial);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                cerr << "ID:" << n1 + ns1 << "," << n2 + ns2 << "," << n3 + ns3
                     << "," << n4 + ns4 << "," << n5 << endl;
                cerr << "ID:Number of isotopes 5. Polynomial size = "
                     << F_Polynomial[k].size() << endl;
            } else if (Element_Composition[k].size() ==
                       6)  // Handles the case of six isotopic elements
            {
                if (Element_Composition[k][0].atoms < n_atoms) {
                    nc_add = 10;
                } else {
                    nc_add = nc_add_value;
                }

                int n1, n2, n3, n4, n5, n6, ns1, ns2, ns3, ns4, ns5, ns6, l, m,
                    n, p;

                // computing mean and standard deviation
                n1 = int(Element_Composition[k][0].atoms *
                         Element_Composition[k][0].prob);
                n2 = int(Element_Composition[k][0].atoms *
                         Element_Composition[k][1].prob);
                n3 = int(Element_Composition[k][0].atoms *
                         Element_Composition[k][2].prob);
                n4 = int(Element_Composition[k][0].atoms *
                         Element_Composition[k][3].prob);
                n5 = int(Element_Composition[k][0].atoms *
                         Element_Composition[k][4].prob);
                n6 = int(Element_Composition[k][0].atoms *
                         Element_Composition[k][5].prob);

                ns1 = int(
                    ceil(nc_add +
                         nc * sqrt(Element_Composition[k][0].atoms *
                                   Element_Composition[k][0].prob *
                                   (1.0 - Element_Composition[k][0].prob))));
                ns2 = int(
                    ceil(nc_add +
                         nc * sqrt(Element_Composition[k][1].atoms *
                                   Element_Composition[k][1].prob *
                                   (1.0 - Element_Composition[k][1].prob))));
                ns3 = int(
                    ceil(nc_add +
                         nc * sqrt(Element_Composition[k][2].atoms *
                                   Element_Composition[k][2].prob *
                                   (1.0 - Element_Composition[k][2].prob))));
                ns4 = int(
                    ceil(nc_add +
                         nc * sqrt(Element_Composition[k][3].atoms *
                                   Element_Composition[k][3].prob *
                                   (1.0 - Element_Composition[k][3].prob))));
                ns5 = int(
                    ceil(nc_add +
                         nc * sqrt(Element_Composition[k][4].atoms *
                                   Element_Composition[k][4].prob *
                                   (1.0 - Element_Composition[k][4].prob))));
                ns6 = int(
                    ceil(nc_add +
                         nc * sqrt(Element_Composition[k][5].atoms *
                                   Element_Composition[k][5].prob *
                                   (1.0 - Element_Composition[k][5].prob))));

                n6 = n6 + ns6;

                n_polynomial_terms = log(n1 + ns1) + log(n2 + ns2) +
                                     log(n3 + ns3) + log(n4 + ns4) +
                                     log(n5 + ns5) + log(n6) + log(pow(2, 6));
                if (n_polynomial_terms >
                    max_polynomial_size)  // checking polynomial size
                {
                    element_composition.push_back(Element_Composition[k]);
                    FT_Fine_Grained_ID(element_composition, T_Polynomial,
                                       FINE_RESOLUTION);
                    for (i = 0; i < int(T_Polynomial.size()); i++) {
                        if (T_Polynomial[i].power > 0) {
                            T_Polynomial[i].power =
                                T_Polynomial[i].power / MW_RESOLUTION;
                            F_Polynomial[k].push_back(T_Polynomial[i]);
                        }
                    }
                    element_composition.clear();
                    T_Polynomial.clear();
                } else {
                    for (i = n1 + ns1; i >= n1 - ns1 && i >= 0; i--) {
                        for (j = n2 + ns2; j >= n2 - ns2 && j >= 0; j--) {
                            for (m = n3 + ns3; m >= n3 - ns3 && m >= 0; m--) {
                                for (n = n4 + ns4; n >= n4 - ns4 && n >= 0;
                                     n--) {
                                    for (p = n5 + ns5; p >= n5 - ns5 && p >= 0;
                                         p--) {
                                        l = Element_Composition[k][0].atoms -
                                            (i + j + m + n + p);
                                        if ((l >= 0) && (l <= n6)) {
                                            prob =
                                                FACTR_LN(
                                                    Element_Composition[k][0]
                                                        .atoms) -
                                                FACTR_LN(i) - FACTR_LN(j) -
                                                FACTR_LN(m) - FACTR_LN(n) -
                                                FACTR_LN(p) - FACTR_LN(l) +
                                                i * (Element_Composition[k][0]
                                                         .log_prob) +
                                                j * (Element_Composition[k][1]
                                                         .log_prob) +
                                                m * (Element_Composition[k][2]
                                                         .log_prob) +
                                                n * (Element_Composition[k][3]
                                                         .log_prob) +
                                                p * (Element_Composition[k][4]
                                                         .log_prob) +
                                                l * (Element_Composition[k][5]
                                                         .log_prob);

                                            prob = exp(prob);

                                            if (prob >= min_prob) {
                                                temp_Polynomial.power =
                                                    i * Element_Composition
                                                            [k][0]
                                                                .power +
                                                    j * Element_Composition
                                                            [k][1]
                                                                .power +
                                                    m * Element_Composition
                                                            [k][2]
                                                                .power +
                                                    n * Element_Composition
                                                            [k][3]
                                                                .power +
                                                    p * Element_Composition
                                                            [k][4]
                                                                .power +
                                                    l * Element_Composition
                                                            [k][5]
                                                                .power;
                                                temp_Polynomial.prob = prob;

                                                F_Polynomial[k].push_back(
                                                    temp_Polynomial);
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                cerr << "ID:" << n1 + ns1 << "," << n2 + ns2 << "," << n3 + ns3
                     << "," << n4 + ns4 << "," << n5 + ns5 << "," << n6 << endl;
                cerr << "ID:Number of isotopes 6. Polynomial size = "
                     << F_Polynomial[k].size() << endl;
            } else if (Element_Composition[k].size() ==
                       7)  // Handles the case of seven isotopic elements
            {
                if (Element_Composition[k][0].atoms < n_atoms) {
                    nc_add = 10;
                } else {
                    nc_add = nc_add_value;
                }

                int n1, n2, n3, n4, n5, n6, n7, ns1, ns2, ns3, ns4, ns5, ns6,
                    ns7, l, m, n, p, q;

                // computing mean and standard deviation
                n1 = int(Element_Composition[k][0].atoms *
                         Element_Composition[k][0].prob);
                n2 = int(Element_Composition[k][0].atoms *
                         Element_Composition[k][1].prob);
                n3 = int(Element_Composition[k][0].atoms *
                         Element_Composition[k][2].prob);
                n4 = int(Element_Composition[k][0].atoms *
                         Element_Composition[k][3].prob);
                n5 = int(Element_Composition[k][0].atoms *
                         Element_Composition[k][4].prob);
                n6 = int(Element_Composition[k][0].atoms *
                         Element_Composition[k][5].prob);
                n7 = int(Element_Composition[k][0].atoms *
                         Element_Composition[k][6].prob);

                ns1 = int(
                    ceil(nc_add +
                         nc * sqrt(Element_Composition[k][0].atoms *
                                   Element_Composition[k][0].prob *
                                   (1.0 - Element_Composition[k][0].prob))));
                ns2 = int(
                    ceil(nc_add +
                         nc * sqrt(Element_Composition[k][1].atoms *
                                   Element_Composition[k][1].prob *
                                   (1.0 - Element_Composition[k][1].prob))));
                ns3 = int(
                    ceil(nc_add +
                         nc * sqrt(Element_Composition[k][2].atoms *
                                   Element_Composition[k][2].prob *
                                   (1.0 - Element_Composition[k][2].prob))));
                ns4 = int(
                    ceil(nc_add +
                         nc * sqrt(Element_Composition[k][3].atoms *
                                   Element_Composition[k][3].prob *
                                   (1.0 - Element_Composition[k][3].prob))));
                ns5 = int(
                    ceil(nc_add +
                         nc * sqrt(Element_Composition[k][4].atoms *
                                   Element_Composition[k][4].prob *
                                   (1.0 - Element_Composition[k][4].prob))));
                ns6 = int(
                    ceil(nc_add +
                         nc * sqrt(Element_Composition[k][5].atoms *
                                   Element_Composition[k][5].prob *
                                   (1.0 - Element_Composition[k][5].prob))));
                ns7 = int(
                    ceil(nc_add +
                         nc * sqrt(Element_Composition[k][6].atoms *
                                   Element_Composition[k][6].prob *
                                   (1.0 - Element_Composition[k][6].prob))));

                n7 = n7 + ns7;
                n_polynomial_terms = log(n1 + ns1) + log(n2 + ns2) +
                                     log(n3 + ns3) + log(n4 + ns4) +
                                     log(n5 + ns5) + log(n6 + ns6) + log(n7) +
                                     log(pow(2, 7));
                if (n_polynomial_terms >
                    max_polynomial_size)  // checking polynomial size
                {
                    element_composition.push_back(Element_Composition[k]);
                    FT_Fine_Grained_ID(element_composition, T_Polynomial,
                                       FINE_RESOLUTION);
                    for (i = 0; i < int(T_Polynomial.size()); i++) {
                        if (T_Polynomial[i].power > 0) {
                            T_Polynomial[i].power =
                                T_Polynomial[i].power / MW_RESOLUTION;
                            F_Polynomial[k].push_back(T_Polynomial[i]);
                        }
                    }
                    element_composition.clear();
                    T_Polynomial.clear();
                } else {
                    for (i = n1 + ns1; i >= n1 - ns1 && i >= 0; i--) {
                        for (j = n2 + ns2; j >= n2 - ns2 && j >= 0; j--) {
                            for (m = n3 + ns3; m >= n3 - ns3 && m >= 0; m--) {
                                for (n = n4 + ns4; n >= n4 - ns4 && n >= 0;
                                     n--) {
                                    for (p = n5 + ns5; p >= n5 - ns5 && p >= 0;
                                         p--) {
                                        for (q = n6 + ns6;
                                             q >= n6 - ns6 && q >= 0; q--) {
                                            l = Element_Composition[k][0]
                                                    .atoms -
                                                (i + j + m + n + p + q);
                                            if ((l >= 0) && (l <= n7)) {
                                                prob =
                                                    FACTR_LN(Element_Composition
                                                                 [k][0]
                                                                     .atoms) -
                                                    FACTR_LN(i) - FACTR_LN(j) -
                                                    FACTR_LN(m) - FACTR_LN(n) -
                                                    FACTR_LN(p) - FACTR_LN(q) -
                                                    FACTR_LN(l) +
                                                    i * (Element_Composition
                                                             [k][0]
                                                                 .log_prob) +
                                                    j * (Element_Composition
                                                             [k][1]
                                                                 .log_prob) +
                                                    m * (Element_Composition
                                                             [k][2]
                                                                 .log_prob) +
                                                    n * (Element_Composition
                                                             [k][3]
                                                                 .log_prob) +
                                                    p * (Element_Composition
                                                             [k][4]
                                                                 .log_prob) +
                                                    q * (Element_Composition
                                                             [k][5]
                                                                 .log_prob) +
                                                    l * (Element_Composition
                                                             [k][6]
                                                                 .log_prob);

                                                prob = exp(prob);
                                                if (prob >= min_prob) {
                                                    temp_Polynomial.power =
                                                        i * Element_Composition
                                                                [k][0]
                                                                    .power +
                                                        j * Element_Composition
                                                                [k][1]
                                                                    .power +
                                                        m * Element_Composition
                                                                [k][2]
                                                                    .power +
                                                        n * Element_Composition
                                                                [k][3]
                                                                    .power +
                                                        p * Element_Composition
                                                                [k][4]
                                                                    .power +
                                                        q * Element_Composition
                                                                [k][5]
                                                                    .power +
                                                        l * Element_Composition
                                                                [k][6]
                                                                    .power;

                                                    temp_Polynomial.prob = prob;

                                                    F_Polynomial[k].push_back(
                                                        temp_Polynomial);
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                cerr << "ID:" << n1 + ns1 << "," << n2 + ns2 << "," << n3 + ns3
                     << "," << n4 + ns4 << "," << n5 + ns5 << "," << n6 + ns6
                     << "," << n7 << endl;
                cerr << "ID:Number of isotopes 7. Polynomial size = "
                     << F_Polynomial[k].size() << endl;
            } else if (Element_Composition[k].size() ==
                       8)  // Handles the case of eight isotopic elements
            {
                if (Element_Composition[k][0].atoms < n_atoms) {
                    nc_add = 10;
                } else {
                    nc_add = nc_add_value;
                }

                int n1, n2, n3, n4, n5, n6, n7, n8, ns1, ns2, ns3, ns4, ns5,
                    ns6, ns7, ns8, l, m, n, p, q, r;

                // computing mean and standard deviation
                n1 = int(Element_Composition[k][0].atoms *
                         Element_Composition[k][0].prob);
                n2 = int(Element_Composition[k][0].atoms *
                         Element_Composition[k][1].prob);
                n3 = int(Element_Composition[k][0].atoms *
                         Element_Composition[k][2].prob);
                n4 = int(Element_Composition[k][0].atoms *
                         Element_Composition[k][3].prob);
                n5 = int(Element_Composition[k][0].atoms *
                         Element_Composition[k][4].prob);
                n6 = int(Element_Composition[k][0].atoms *
                         Element_Composition[k][5].prob);
                n7 = int(Element_Composition[k][0].atoms *
                         Element_Composition[k][6].prob);
                n8 = int(Element_Composition[k][0].atoms *
                         Element_Composition[k][7].prob);

                ns1 = int(
                    ceil(nc_add +
                         nc * sqrt(Element_Composition[k][0].atoms *
                                   Element_Composition[k][0].prob *
                                   (1.0 - Element_Composition[k][0].prob))));
                ns2 = int(
                    ceil(nc_add +
                         nc * sqrt(Element_Composition[k][1].atoms *
                                   Element_Composition[k][1].prob *
                                   (1.0 - Element_Composition[k][1].prob))));
                ns3 = int(
                    ceil(nc_add +
                         nc * sqrt(Element_Composition[k][2].atoms *
                                   Element_Composition[k][2].prob *
                                   (1.0 - Element_Composition[k][2].prob))));
                ns4 = int(
                    ceil(nc_add +
                         nc * sqrt(Element_Composition[k][3].atoms *
                                   Element_Composition[k][3].prob *
                                   (1.0 - Element_Composition[k][3].prob))));
                ns5 = int(
                    ceil(nc_add +
                         nc * sqrt(Element_Composition[k][4].atoms *
                                   Element_Composition[k][4].prob *
                                   (1.0 - Element_Composition[k][4].prob))));
                ns6 = int(
                    ceil(nc_add +
                         nc * sqrt(Element_Composition[k][5].atoms *
                                   Element_Composition[k][5].prob *
                                   (1.0 - Element_Composition[k][5].prob))));
                ns7 = int(
                    ceil(nc_add +
                         nc * sqrt(Element_Composition[k][6].atoms *
                                   Element_Composition[k][6].prob *
                                   (1.0 - Element_Composition[k][6].prob))));
                ns8 = int(
                    ceil(nc_add +
                         nc * sqrt(Element_Composition[k][7].atoms *
                                   Element_Composition[k][7].prob *
                                   (1.0 - Element_Composition[k][7].prob))));
                n8 = n8 + ns8;

                n_polynomial_terms = log(n1 + ns1) + log(n2 + ns2) +
                                     log(n3 + ns3) + log(n4 + ns4) +
                                     log(n5 + ns5) + log(n6 + ns6) +
                                     log(n7 + ns7) + log(n8) + log(pow(2, 8));
                if (n_polynomial_terms >
                    max_polynomial_size)  // checking polynomial size
                {
                    element_composition.push_back(Element_Composition[k]);
                    FT_Fine_Grained_ID(element_composition, T_Polynomial,
                                       FINE_RESOLUTION);
                    for (i = 0; i < int(T_Polynomial.size()); i++) {
                        if (T_Polynomial[i].power > 0) {
                            T_Polynomial[i].power =
                                T_Polynomial[i].power / MW_RESOLUTION;
                            F_Polynomial[k].push_back(T_Polynomial[i]);
                        }
                    }
                    element_composition.clear();
                    T_Polynomial.clear();
                } else {
                    for (i = n1 + ns1; i >= n1 - ns1 && i >= 0; i--) {
                        for (j = n2 + ns2; j >= n2 - ns2 && j >= 0; j--) {
                            for (m = n3 + ns3; m >= n3 - ns3 && m >= 0; m--) {
                                for (n = n4 + ns4; n >= n4 - ns4 && n >= 0;
                                     n--) {
                                    for (p = n5 + ns5; p >= n5 - ns5 && p >= 0;
                                         p--) {
                                        for (q = n6 + ns6;
                                             q >= n6 - ns6 && q >= 0; q--) {
                                            for (r = n7 + ns7;
                                                 r >= n7 - ns7 && r >= 0; r--) {
                                                l = Element_Composition[k][0]
                                                        .atoms -
                                                    (i + j + m + n + p + q + r);
                                                if ((l >= 0) && (l <= n8)) {
                                                    prob =
                                                        FACTR_LN(
                                                            Element_Composition
                                                                [k][0]
                                                                    .atoms) -
                                                        FACTR_LN(i) -
                                                        FACTR_LN(j) -
                                                        FACTR_LN(m) -
                                                        FACTR_LN(n) -
                                                        FACTR_LN(p) -
                                                        FACTR_LN(q) -
                                                        FACTR_LN(r) -
                                                        FACTR_LN(l) +
                                                        i * (Element_Composition
                                                                 [k][0]
                                                                     .log_prob) +
                                                        j * (Element_Composition
                                                                 [k][1]
                                                                     .log_prob) +
                                                        m * (Element_Composition
                                                                 [k][2]
                                                                     .log_prob) +
                                                        n * (Element_Composition
                                                                 [k][3]
                                                                     .log_prob) +
                                                        p * (Element_Composition
                                                                 [k][4]
                                                                     .log_prob) +
                                                        q * (Element_Composition
                                                                 [k][5]
                                                                     .log_prob) +
                                                        r * (Element_Composition
                                                                 [k][6]
                                                                     .log_prob) +
                                                        l * (Element_Composition
                                                                 [k][7]
                                                                     .log_prob);

                                                    prob = exp(prob);

                                                    if (prob >= min_prob) {
                                                        temp_Polynomial.power =
                                                            i * Element_Composition
                                                                    [k][0]
                                                                        .power +
                                                            j * Element_Composition
                                                                    [k][1]
                                                                        .power +
                                                            m * Element_Composition
                                                                    [k][2]
                                                                        .power +
                                                            n * Element_Composition
                                                                    [k][3]
                                                                        .power +
                                                            p * Element_Composition
                                                                    [k][4]
                                                                        .power +
                                                            q * Element_Composition
                                                                    [k][5]
                                                                        .power +
                                                            r * Element_Composition
                                                                    [k][6]
                                                                        .power +
                                                            l * Element_Composition
                                                                    [k][7]
                                                                        .power;

                                                        temp_Polynomial.prob =
                                                            prob;

                                                        F_Polynomial[k]
                                                            .push_back(
                                                                temp_Polynomial);
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                cerr << "ID:" << n1 + ns1 << "," << n2 + ns2 << "," << n3 + ns3
                     << "," << n4 + ns4 << "," << n5 + ns5 << "," << n6 + ns6
                     << "," << n7 + ns7 << "," << n8 << endl;
                cerr << "ID:Number of isotopes 8. Polynomial size = "
                     << F_Polynomial[k].size() << endl;
            } else if (Element_Composition[k].size() ==
                       9)  // Handles the case of nine isotopic elements
            {
                if (Element_Composition[k][0].atoms < n_atoms) {
                    nc_add = 10;
                } else {
                    nc_add = nc_add_value;
                }

                int n1, n2, n3, n4, n5, n6, n7, n8, n9, ns1, ns2, ns3, ns4, ns5,
                    ns6, ns7, ns8, ns9, l, m, n, p, q, r, s;

                // computing mean and standard deviation
                n1 = int(Element_Composition[k][0].atoms *
                         Element_Composition[k][0].prob);
                n2 = int(Element_Composition[k][0].atoms *
                         Element_Composition[k][1].prob);
                n3 = int(Element_Composition[k][0].atoms *
                         Element_Composition[k][2].prob);
                n4 = int(Element_Composition[k][0].atoms *
                         Element_Composition[k][3].prob);
                n5 = int(Element_Composition[k][0].atoms *
                         Element_Composition[k][4].prob);
                n6 = int(Element_Composition[k][0].atoms *
                         Element_Composition[k][5].prob);
                n7 = int(Element_Composition[k][0].atoms *
                         Element_Composition[k][6].prob);
                n8 = int(Element_Composition[k][0].atoms *
                         Element_Composition[k][7].prob);
                n9 = int(Element_Composition[k][0].atoms *
                         Element_Composition[k][8].prob);

                ns1 = int(
                    ceil(nc_add +
                         nc * sqrt(Element_Composition[k][0].atoms *
                                   Element_Composition[k][0].prob *
                                   (1.0 - Element_Composition[k][0].prob))));
                ns2 = int(
                    ceil(nc_add +
                         nc * sqrt(Element_Composition[k][1].atoms *
                                   Element_Composition[k][1].prob *
                                   (1.0 - Element_Composition[k][1].prob))));
                ns3 = int(
                    ceil(nc_add +
                         nc * sqrt(Element_Composition[k][2].atoms *
                                   Element_Composition[k][2].prob *
                                   (1.0 - Element_Composition[k][2].prob))));
                ns4 = int(
                    ceil(nc_add +
                         nc * sqrt(Element_Composition[k][3].atoms *
                                   Element_Composition[k][3].prob *
                                   (1.0 - Element_Composition[k][3].prob))));
                ns5 = int(
                    ceil(nc_add +
                         nc * sqrt(Element_Composition[k][4].atoms *
                                   Element_Composition[k][4].prob *
                                   (1.0 - Element_Composition[k][4].prob))));
                ns6 = int(
                    ceil(nc_add +
                         nc * sqrt(Element_Composition[k][5].atoms *
                                   Element_Composition[k][5].prob *
                                   (1.0 - Element_Composition[k][5].prob))));
                ns7 = int(
                    ceil(nc_add +
                         nc * sqrt(Element_Composition[k][6].atoms *
                                   Element_Composition[k][6].prob *
                                   (1.0 - Element_Composition[k][6].prob))));
                ns8 = int(
                    ceil(nc_add +
                         nc * sqrt(Element_Composition[k][7].atoms *
                                   Element_Composition[k][7].prob *
                                   (1.0 - Element_Composition[k][7].prob))));
                ns9 = int(
                    ceil(nc_add +
                         nc * sqrt(Element_Composition[k][8].atoms *
                                   Element_Composition[k][8].prob *
                                   (1.0 - Element_Composition[k][8].prob))));

                n9 = n9 + ns9;
                n_polynomial_terms =
                    log(n1 + ns1) + log(n2 + ns2) + log(n3 + ns3) +
                    log(n4 + ns4) + log(n5 + ns5) + log(n6 + ns6) +
                    log(n7 + ns7) + log(n8 + ns8) + log(n9) + log(pow(2, 9));
                if (n_polynomial_terms >
                    max_polynomial_size)  // checking polynomial size
                {
                    element_composition.push_back(Element_Composition[k]);
                    FT_Fine_Grained_ID(element_composition, T_Polynomial,
                                       FINE_RESOLUTION);
                    for (i = 0; i < int(T_Polynomial.size()); i++) {
                        if (T_Polynomial[i].power > 0) {
                            T_Polynomial[i].power =
                                T_Polynomial[i].power / MW_RESOLUTION;
                            F_Polynomial[k].push_back(T_Polynomial[i]);
                        }
                    }
                    element_composition.clear();
                    T_Polynomial.clear();
                } else {
                    for (i = n1 + ns1; i >= n1 - ns1 && i >= 0; i--) {
                        for (j = n2 + ns2; j >= n2 - ns2 && j >= 0; j--) {
                            for (m = n3 + ns3; m >= n3 - ns3 && m >= 0; m--) {
                                for (n = n4 + ns4; n >= n4 - ns4 && n >= 0;
                                     n--) {
                                    for (p = n5 + ns5; p >= n5 - ns5 && p >= 0;
                                         p--) {
                                        for (q = n6 + ns6;
                                             q >= n6 - ns6 && q >= 0; q--) {
                                            for (r = n7 + ns7;
                                                 r >= n7 - ns7 && r >= 0; r--) {
                                                for (s = n8 + ns8;
                                                     s >= n8 - ns8 && s >= 0;
                                                     s--) {
                                                    l = Element_Composition
                                                            [k][0]
                                                                .atoms -
                                                        (i + j + m + n + p + q +
                                                         r + s);
                                                    if ((l >= 0) && (l <= n9)) {
                                                        prob =
                                                            FACTR_LN(
                                                                Element_Composition
                                                                    [k][0]
                                                                        .atoms) -
                                                            FACTR_LN(i) -
                                                            FACTR_LN(j) -
                                                            FACTR_LN(m) -
                                                            FACTR_LN(n) -
                                                            FACTR_LN(p) -
                                                            FACTR_LN(q) -
                                                            FACTR_LN(r) -
                                                            FACTR_LN(s) -
                                                            FACTR_LN(l) +
                                                            i * (Element_Composition
                                                                     [k][0]
                                                                         .log_prob) +
                                                            j * (Element_Composition
                                                                     [k][1]
                                                                         .log_prob) +
                                                            m * (Element_Composition
                                                                     [k][2]
                                                                         .log_prob) +
                                                            n * (Element_Composition
                                                                     [k][3]
                                                                         .log_prob) +
                                                            p * (Element_Composition
                                                                     [k][4]
                                                                         .log_prob) +
                                                            q * (Element_Composition
                                                                     [k][5]
                                                                         .log_prob) +
                                                            r * (Element_Composition
                                                                     [k][6]
                                                                         .log_prob) +
                                                            s * (Element_Composition
                                                                     [k][7]
                                                                         .log_prob) +
                                                            l * (Element_Composition
                                                                     [k][8]
                                                                         .log_prob);

                                                        prob = exp(prob);

                                                        if (prob >= min_prob) {
                                                            temp_Polynomial
                                                                .power =
                                                                i * Element_Composition
                                                                        [k][0]
                                                                            .power +
                                                                j * Element_Composition
                                                                        [k][1]
                                                                            .power +
                                                                m * Element_Composition
                                                                        [k][2]
                                                                            .power +
                                                                n * Element_Composition
                                                                        [k][3]
                                                                            .power +
                                                                p * Element_Composition
                                                                        [k][4]
                                                                            .power +
                                                                q * Element_Composition
                                                                        [k][5]
                                                                            .power +
                                                                r * Element_Composition
                                                                        [k][6]
                                                                            .power +
                                                                s * Element_Composition
                                                                        [k][7]
                                                                            .power +
                                                                l * Element_Composition
                                                                        [k][8]
                                                                            .power;

                                                            temp_Polynomial
                                                                .prob = prob;

                                                            F_Polynomial[k].push_back(
                                                                temp_Polynomial);
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                cerr << "ID:" << n1 + ns1 << "," << n2 + ns2 << "," << n3 + ns3
                     << "," << n4 + ns4 << "," << n5 + ns5 << "," << n6 + ns6
                     << "," << n7 + ns7 << "," << n8 + ns8 << "," << n9 << endl;
                cerr << "ID:Number of isotopes 9. Polynomial size = "
                     << F_Polynomial[k].size() << endl;
            } else if (Element_Composition[k].size() ==
                       10)  // Handles the case of ten isotopic elements
            {
                if (Element_Composition[k][0].atoms < n_atoms) {
                    nc_add = 10;
                } else {
                    nc_add = nc_add_value;
                }

                int n1, n2, n3, n4, n5, n6, n7, n8, n9, n10, ns1, ns2, ns3, ns4,
                    ns5, ns6, ns7, ns8, ns9, ns10, l, m, n, p, q, r, s, t;

                // computing mean and standard deviation
                n1 = int(Element_Composition[k][0].atoms *
                         Element_Composition[k][0].prob);
                n2 = int(Element_Composition[k][0].atoms *
                         Element_Composition[k][1].prob);
                n3 = int(Element_Composition[k][0].atoms *
                         Element_Composition[k][2].prob);
                n4 = int(Element_Composition[k][0].atoms *
                         Element_Composition[k][3].prob);
                n5 = int(Element_Composition[k][0].atoms *
                         Element_Composition[k][4].prob);
                n6 = int(Element_Composition[k][0].atoms *
                         Element_Composition[k][5].prob);
                n7 = int(Element_Composition[k][0].atoms *
                         Element_Composition[k][6].prob);
                n8 = int(Element_Composition[k][0].atoms *
                         Element_Composition[k][7].prob);
                n9 = int(Element_Composition[k][0].atoms *
                         Element_Composition[k][8].prob);
                n10 = int(Element_Composition[k][0].atoms *
                          Element_Composition[k][9].prob);

                ns1 = int(
                    ceil(nc_add +
                         nc * sqrt(Element_Composition[k][0].atoms *
                                   Element_Composition[k][0].prob *
                                   (1.0 - Element_Composition[k][0].prob))));
                ns2 = int(
                    ceil(nc_add +
                         nc * sqrt(Element_Composition[k][1].atoms *
                                   Element_Composition[k][1].prob *
                                   (1.0 - Element_Composition[k][1].prob))));
                ns3 = int(
                    ceil(nc_add +
                         nc * sqrt(Element_Composition[k][2].atoms *
                                   Element_Composition[k][2].prob *
                                   (1.0 - Element_Composition[k][2].prob))));
                ns4 = int(
                    ceil(nc_add +
                         nc * sqrt(Element_Composition[k][3].atoms *
                                   Element_Composition[k][3].prob *
                                   (1.0 - Element_Composition[k][3].prob))));
                ns5 = int(
                    ceil(nc_add +
                         nc * sqrt(Element_Composition[k][4].atoms *
                                   Element_Composition[k][4].prob *
                                   (1.0 - Element_Composition[k][4].prob))));
                ns6 = int(
                    ceil(nc_add +
                         nc * sqrt(Element_Composition[k][5].atoms *
                                   Element_Composition[k][5].prob *
                                   (1.0 - Element_Composition[k][5].prob))));
                ns7 = int(
                    ceil(nc_add +
                         nc * sqrt(Element_Composition[k][6].atoms *
                                   Element_Composition[k][6].prob *
                                   (1.0 - Element_Composition[k][6].prob))));
                ns8 = int(
                    ceil(nc_add +
                         nc * sqrt(Element_Composition[k][7].atoms *
                                   Element_Composition[k][7].prob *
                                   (1.0 - Element_Composition[k][7].prob))));
                ns9 = int(
                    ceil(nc_add +
                         nc * sqrt(Element_Composition[k][8].atoms *
                                   Element_Composition[k][8].prob *
                                   (1.0 - Element_Composition[k][8].prob))));
                ns10 = int(
                    ceil(nc_add +
                         nc * sqrt(Element_Composition[k][9].atoms *
                                   Element_Composition[k][9].prob *
                                   (1.0 - Element_Composition[k][9].prob))));

                n10 = n10 + ns10;

                n_polynomial_terms = log(n1 + ns1) + log(n2 + ns2) +
                                     log(n3 + ns3) + log(n4 + ns4) +
                                     log(n5 + ns5) + log(n6 + ns6) +
                                     log(n7 + ns7) + log(n8 + ns8) +
                                     log(n9 + ns9) + log(n10) + log(pow(2, 10));
                if (n_polynomial_terms >
                    max_polynomial_size)  // checking polynomial size
                {
                    element_composition.push_back(Element_Composition[k]);
                    FT_Fine_Grained_ID(element_composition, T_Polynomial,
                                       FINE_RESOLUTION);
                    for (i = 0; i < int(T_Polynomial.size()); i++) {
                        if (T_Polynomial[i].power > 0) {
                            T_Polynomial[i].power =
                                T_Polynomial[i].power / MW_RESOLUTION;
                            F_Polynomial[k].push_back(T_Polynomial[i]);
                        }
                    }
                    element_composition.clear();
                    T_Polynomial.clear();
                } else {
                    for (i = n1 + ns1; i >= n1 - ns1 && i >= 0; i--) {
                        for (j = n2 + ns2; j >= n2 - ns2 && j >= 0; j--) {
                            for (m = n3 + ns3; m >= n3 - ns3 && m >= 0; m--) {
                                for (n = n4 + ns4; n >= n4 - ns4 && n >= 0;
                                     n--) {
                                    for (p = n5 + ns5; p >= n5 - ns5 && p >= 0;
                                         p--) {
                                        for (q = n6 + ns6;
                                             q >= n6 - ns6 && q >= 0; q--) {
                                            for (r = n7 + ns7;
                                                 r >= n7 - ns7 && r >= 0; r--) {
                                                for (s = n8 + ns8;
                                                     s >= n8 - ns8 && s >= 0;
                                                     s--) {
                                                    for (t = n9 + ns9;
                                                         t >= n9 - ns9 &&
                                                         t >= 0;
                                                         t--) {
                                                        l = Element_Composition
                                                                [k][0]
                                                                    .atoms -
                                                            (i + j + m + n + p +
                                                             q + r + s + t);
                                                        if ((l >= 0) &&
                                                            (l <= n10)) {
                                                            prob =
                                                                FACTR_LN(
                                                                    Element_Composition
                                                                        [k][0]
                                                                            .atoms) -
                                                                FACTR_LN(i) -
                                                                FACTR_LN(j) -
                                                                FACTR_LN(m) -
                                                                FACTR_LN(n) -
                                                                FACTR_LN(p) -
                                                                FACTR_LN(q) -
                                                                FACTR_LN(r) -
                                                                FACTR_LN(s) -
                                                                FACTR_LN(t) -
                                                                FACTR_LN(l) +
                                                                i * (Element_Composition
                                                                         [k][0]
                                                                             .log_prob) +
                                                                j * (Element_Composition
                                                                         [k][1]
                                                                             .log_prob) +
                                                                m * (Element_Composition
                                                                         [k][2]
                                                                             .log_prob) +
                                                                n * (Element_Composition
                                                                         [k][3]
                                                                             .log_prob) +
                                                                p * (Element_Composition
                                                                         [k][4]
                                                                             .log_prob) +
                                                                q * (Element_Composition
                                                                         [k][5]
                                                                             .log_prob) +
                                                                r * (Element_Composition
                                                                         [k][6]
                                                                             .log_prob) +
                                                                s * (Element_Composition
                                                                         [k][7]
                                                                             .log_prob) +
                                                                t * (Element_Composition
                                                                         [k][8]
                                                                             .log_prob) +
                                                                l * (Element_Composition
                                                                         [k][9]
                                                                             .log_prob);

                                                            prob = exp(prob);

                                                            if (prob >=
                                                                min_prob) {
                                                                temp_Polynomial
                                                                    .power =
                                                                    i * Element_Composition
                                                                            [k]
                                                                            [0]
                                                                                .power +
                                                                    j * Element_Composition
                                                                            [k]
                                                                            [1]
                                                                                .power +
                                                                    m * Element_Composition
                                                                            [k]
                                                                            [2]
                                                                                .power +
                                                                    n * Element_Composition
                                                                            [k]
                                                                            [3]
                                                                                .power +
                                                                    p * Element_Composition
                                                                            [k]
                                                                            [4]
                                                                                .power +
                                                                    q * Element_Composition
                                                                            [k]
                                                                            [5]
                                                                                .power +
                                                                    r * Element_Composition
                                                                            [k]
                                                                            [6]
                                                                                .power +
                                                                    s * Element_Composition
                                                                            [k]
                                                                            [7]
                                                                                .power +
                                                                    t * Element_Composition
                                                                            [k]
                                                                            [8]
                                                                                .power +
                                                                    l * Element_Composition
                                                                            [k]
                                                                            [9]
                                                                                .power;

                                                                temp_Polynomial
                                                                    .prob =
                                                                    prob;

                                                                F_Polynomial[k]
                                                                    .push_back(
                                                                        temp_Polynomial);
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                cerr << "ID:" << n1 + ns1 << "," << n2 + ns2 << "," << n3 + ns3
                     << "," << n4 + ns4 << "," << n5 + ns5 << "," << n6 + ns6
                     << "," << n7 + ns7 << "," << n8 + ns8 << "," << n9 + ns9
                     << "," << n10 << endl;
                cerr << "ID:Number of isotopes 10. Polynomial size = "
                     << F_Polynomial[k].size() << endl;
            } else  // Handling case of an element having more than 10 Isotopes
            {
                element_composition.push_back(Element_Composition[k]);
                FT_Fine_Grained_ID(element_composition, T_Polynomial,
                                   FINE_RESOLUTION);
                for (i = 0; i < int(T_Polynomial.size()); i++) {
                    if (T_Polynomial[i].power > 0) {
                        T_Polynomial[i].power =
                            T_Polynomial[i].power / MW_RESOLUTION;
                        F_Polynomial[k].push_back(T_Polynomial[i]);
                    }
                }
                element_composition.clear();
                T_Polynomial.clear();
                cerr << "ID:Number of isotopes" << Element_Composition[k].size()
                     << ". Polynomial size = " << F_Polynomial[k].size()
                     << endl;
            }
        }

        cerr << "ID:-----------------------------------------------------------"
                "-----------------------------"
             << endl;
        cerr << "ID:Starting computing cross element polynomial products"
             << endl;
        cerr << "ID:-----------------------------------------------------------"
                "-----------------------------"
             << endl;
        if (k > 1) {
            T_Polynomial = F_Polynomial[0];

            for (k = 1; k < n; k++) {
                cerr << "ID:Multiplication step = " << k << endl;
                cerr << "ID:Size of polynomial multiplying = "
                     << T_Polynomial.size() << " X " << F_Polynomial[k].size()
                     << endl;
                multiplication_fine_final_operation(T_Polynomial,
                                                    F_Polynomial[k]);
            }
        } else {
            T_Polynomial = F_Polynomial[0];
        }
    }
    cerr << "ID:---------------------------------------------------------------"
            "-------------------------"
         << endl;
    return T_Polynomial;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*Multiplying the polynomials from different polynomials for fine isotopic
 * distribution*/
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void MIDAs::multiplication_fine_final_operation(
    vector<struct Polynomial> &T_Polynomial, vector<Polynomial> &F_Polynomial) {
    unsigned int i, j, index, max_index;
    double min_prob;
    double delta_mass;
    double max_mass, min_mass;
    struct Polynomial temp_Polynomial;

    delta_mass = (FINE_RESOLUTION / MW_RESOLUTION);

    min_prob = FINE_MIN_PROB;
    // computing max_index
    i = T_Polynomial.size();
    j = F_Polynomial.size();
    min_mass = T_Polynomial[0].power + F_Polynomial[0].power;
    max_mass = T_Polynomial[i - 1].power + F_Polynomial[j - 1].power;
    max_index = unsigned(fabs(max_mass - min_mass) / delta_mass + 0.5);
    if (max_index > VECTOR_FGID_SIZE)  // resize vector
    {
        j = max_index - VECTOR_FGID_SIZE;
        for (i = 0; i <= j; i++) {
            temp_Polynomial.prob = temp_Polynomial.power = 0;
            FGID_Polynomial.push_back(temp_Polynomial);
        }
        VECTOR_FGID_SIZE = FGID_Polynomial.size();
    }

    for (i = 0; i < T_Polynomial.size(); i++) {
        for (j = 0; j < F_Polynomial.size(); j++) {
            temp_Polynomial.prob = T_Polynomial[i].prob * F_Polynomial[j].prob;

            if (temp_Polynomial.prob > min_prob) {
                temp_Polynomial.power =
                    T_Polynomial[i].power + F_Polynomial[j].power;
                index = unsigned(
                    fabs(temp_Polynomial.power - min_mass) / delta_mass + 0.5);
                FGID_Polynomial[index].prob =
                    FGID_Polynomial[index].prob + temp_Polynomial.prob;
                FGID_Polynomial[index].power =
                    FGID_Polynomial[index].power +
                    temp_Polynomial.prob * temp_Polynomial.power;
            }
        }
    }

    index = T_Polynomial.size();
    j = 0;
    for (i = 0; i < FGID_Polynomial.size(); i++) {
        if (FGID_Polynomial[i].prob != 0) {
            if (j < index) {
                T_Polynomial[j].prob = FGID_Polynomial[i].prob;
                T_Polynomial[j].power =
                    FGID_Polynomial[i].power / FGID_Polynomial[i].prob;
                j = j + 1;
            } else {
                temp_Polynomial.prob = FGID_Polynomial[i].prob;
                temp_Polynomial.power =
                    FGID_Polynomial[i].power / FGID_Polynomial[i].prob;
                T_Polynomial.push_back(temp_Polynomial);
            }
        }

        FGID_Polynomial[i].prob = 0;
        FGID_Polynomial[i].power = 0;
    }

    if (j < index) {
        T_Polynomial.erase(T_Polynomial.begin() + j, T_Polynomial.end());
    }
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*Computing Coarse-Grained Isotopic Distribution using polynomial
 * multiplication.*/
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
vector<struct Polynomial> MIDAs::multiply_polynomial(
    vector<vector<struct Composition> > &Element_Composition) {
    int i, k, n;
    int N, M, R;
    double min_prob = 1.0;
    struct Polynomial temp_Polynomial;
    vector<struct Polynomial> T_Polynomial;
    vector<struct Isotopic_Distribution> M_Polynomial;
    vector<vector<struct Composition> > element_composition;
    cout.precision(12);

    // Count number of unique elements
    n = 0;
    for (k = 0; k < int(Element_Composition.size()); k++) {
        if (Element_Composition[k].size() > 0) {
            n = n + 1;
        }
    }

    // Start individual polynomial multiplication for individual elements
    vector<vector<struct Polynomial> > F_Polynomial(n);
    cerr << "ID:Polynomial expansion in progress" << endl;
    cerr << "ID:---------------------------------------------------------------"
            "-------------------------"
         << endl;
    N = 0;
    R = 10;
    for (k = 0; k < n; k++) {
        N = 0;
        R = 10;
        T_Polynomial.clear();
        temp_Polynomial.power = 0;
        temp_Polynomial.prob = 1.0;
        T_Polynomial.push_back(temp_Polynomial);

        cerr << "ID:"
             << "Element = " << Element_Composition[k][0].element
             << ",number of atoms = " << Element_Composition[k][0].atoms
             << ",number of isotopics = " << Element_Composition[k].size()
             << endl;

        if (Element_Composition[k].size() >
            10)  // handling case of more than 10 isotopic elements
        {
            element_composition.push_back(Element_Composition[k]);
            T_Polynomial.clear();
            FT_Coarse_Grained_ID(element_composition, M_Polynomial,
                                 COARSE_RESOLUTION);
            for (i = 0; i < int(M_Polynomial.size()); i++) {
                if (M_Polynomial[i].mw > 0) {
                    temp_Polynomial.power = M_Polynomial[i].mw / MW_RESOLUTION;
                    temp_Polynomial.prob = M_Polynomial[i].prob;
                    if (temp_Polynomial.prob < min_prob) {
                        min_prob = temp_Polynomial.prob;
                    }
                    T_Polynomial.push_back(temp_Polynomial);
                }
            }
            element_composition.clear();
            M_Polynomial.clear();

        } else {
            N = int(1.0 * int(Element_Composition[k][0].atoms) / (1.0 * R));
            M = int(fmod(1.0 * int(Element_Composition[k][0].atoms), 1.0 * R));
            cerr << "ID:" << Element_Composition[k][0].atoms << " mod " << R
                 << " = " << M << ",quotient = " << N << endl;

            if (N > 0) {
                cerr << "ID:Start quotient multiplication, vector size = "
                     << T_Polynomial.size() << endl;
                for (i = 0; i < N; i++)  // multiplying  power of N polynomial
                {
                    multiplication_operation(Element_Composition[k],
                                             T_Polynomial);
                }
                F_Polynomial[k] = T_Polynomial;

                cerr << "ID:Finish quotient multiplication" << endl;
                cerr << "ID:Start divisor multiplication, vector size = "
                     << T_Polynomial.size() << endl;

                for (i = 1; i < R; i++)  // multiplying power of R
                {
                    multiplication_final_operation(F_Polynomial[k],
                                                   T_Polynomial);
                }
                cerr << "ID:Finish divisor multiplication" << endl;
            }
            if (M > 0)  // multiplying remainder M of element composition to
                        // polynomial
            {
                cerr << "ID:Start remainder multiplication, vector size =  "
                     << T_Polynomial.size() << endl;
                for (i = 0; i < M; i++) {
                    multiplication_operation(Element_Composition[k],
                                             T_Polynomial);
                }
                cerr << "ID:Finish remainder multiplication" << endl;
            }
        }

        F_Polynomial[k] = T_Polynomial;
        cerr << "ID:Finish computing polynomial for element = "
             << Element_Composition[k][0].element
             << ", polynomial size = " << F_Polynomial[k].size() << endl;
        cerr << "ID:-----------------------------------------------------------"
                "-----------------------------"
             << endl;
    }

    cerr << "ID:Start computing cross element polynomial products" << endl;
    if (min_prob < 1.0) {
        if (min_prob > MIN_PROB) {
            MIN_PROB = min_prob;
        }
    }  // In case FT is invoked
    cerr << "ID:Minimum probability = " << MIN_PROB << endl;
    if (k > 1) {
        cerr << "ID:Multiplication" << endl;
        T_Polynomial.clear();
        T_Polynomial = F_Polynomial[0];

        for (k = 1; k < n; k++) {
            cerr << "ID:Multiplication step = " << k << endl;
            cerr << "ID:Size of polynomial multiplying = "
                 << T_Polynomial.size() << " X " << F_Polynomial[k].size()
                 << endl;
            multiplication_final_operation(F_Polynomial[k], T_Polynomial);
        }
    }
    cerr << "ID:Finish computing cross element polynomial products" << endl;
    cerr << "ID:---------------------------------------------------------------"
            "-------------------------"
         << endl;

    return T_Polynomial;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*Multiplying polynomials from different elements, i.e., C.N,H,O,S,P and etc*/
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void MIDAs::multiplication_final_operation(
    vector<struct Polynomial> &E_Polynomial,
    vector<struct Polynomial> &T_Polynomial) {
    int flag_size;
    unsigned int i, j, index, max_index;
    double delta_mass;
    double max_mass, min_mass;
    struct Polynomial temp_Polynomial;

    delta_mass = (COARSE_RESOLUTION / MW_RESOLUTION);
    // computing max_index
    i = T_Polynomial.size();
    j = E_Polynomial.size();
    min_mass = T_Polynomial[0].power + E_Polynomial[0].power;
    max_mass = T_Polynomial[i - 1].power + E_Polynomial[j - 1].power;
    max_index = unsigned(fabs(max_mass - min_mass) / delta_mass + 0.5);

    if (max_index >= VECTOR_CGID_SIZE)  // resize vector
    {
        j = max_index - VECTOR_CGID_SIZE;
        for (i = 0; i <= j; i++) {
            temp_Polynomial.prob = temp_Polynomial.power = 0;
            CGID_Polynomial.push_back(temp_Polynomial);
        }
        VECTOR_CGID_SIZE = CGID_Polynomial.size();
    }

    flag_size = int(E_Polynomial.size() * T_Polynomial.size());
    index = 0;
    for (i = 0; i < E_Polynomial.size(); i++) {
        for (j = 0; j < T_Polynomial.size(); j++) {
            temp_Polynomial.prob = E_Polynomial[i].prob * T_Polynomial[j].prob;
            if (temp_Polynomial.prob > MIN_PROB) {
                temp_Polynomial.power =
                    E_Polynomial[i].power + T_Polynomial[j].power;

                if (flag_size > FINE_GRID) {
                    index = unsigned(fabs(temp_Polynomial.power - min_mass) /
                                         delta_mass +
                                     0.5);
                    CGID_Polynomial[index].prob =
                        CGID_Polynomial[index].prob + temp_Polynomial.prob;
                    CGID_Polynomial[index].power =
                        CGID_Polynomial[index].power +
                        temp_Polynomial.prob * temp_Polynomial.power;
                } else {
                    CGID_Polynomial[index].prob =
                        CGID_Polynomial[index].prob + temp_Polynomial.prob;
                    CGID_Polynomial[index].power =
                        CGID_Polynomial[index].power +
                        temp_Polynomial.prob * temp_Polynomial.power;
                    index = index + 1;
                }
            }
        }
    }

    index = T_Polynomial.size();
    j = 0;
    for (i = 0; i < CGID_Polynomial.size(); i++) {
        if (CGID_Polynomial[i].prob != 0) {
            if (j < index) {
                T_Polynomial[j].prob = CGID_Polynomial[i].prob;
                T_Polynomial[j].power =
                    CGID_Polynomial[i].power / CGID_Polynomial[i].prob;
                j = j + 1;
            } else {
                temp_Polynomial.prob = CGID_Polynomial[i].prob;
                temp_Polynomial.power =
                    CGID_Polynomial[i].power / CGID_Polynomial[i].prob;
                T_Polynomial.push_back(temp_Polynomial);
            }
        }

        CGID_Polynomial[i].prob = 0;
        CGID_Polynomial[i].power = 0;
    }
    if (j < index) {
        T_Polynomial.erase(T_Polynomial.begin() + j, T_Polynomial.end());
    }
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*Computing an element's polynomial or euivalent to the element isotopic
 * distirbution.*/
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void MIDAs::multiplication_operation(vector<struct Composition> &E_Polynomial,
                                     vector<struct Polynomial> &T_Polynomial) {
    int flag_size;
    unsigned int i, j, index, max_index;
    double delta_mass;
    double max_mass, min_mass;
    struct Polynomial temp_Polynomial;

    delta_mass = (COARSE_RESOLUTION / MW_RESOLUTION);
    // computing max_index
    i = T_Polynomial.size();
    j = E_Polynomial.size();
    min_mass = T_Polynomial[0].power + E_Polynomial[0].power;
    max_mass = T_Polynomial[i - 1].power + E_Polynomial[j - 1].power;
    max_index = unsigned(fabs(max_mass - min_mass) / delta_mass + 0.5);

    if (max_index >= VECTOR_CGID_SIZE)  // resize vector
    {
        j = max_index - VECTOR_CGID_SIZE;
        for (i = 0; i <= j; i++) {
            temp_Polynomial.prob = temp_Polynomial.power = 0;
            CGID_Polynomial.push_back(temp_Polynomial);
        }
        VECTOR_CGID_SIZE = CGID_Polynomial.size();
    }

    flag_size = int(E_Polynomial.size() * T_Polynomial.size());
    index = 0;
    for (i = 0; i < E_Polynomial.size(); i++) {
        // cout<<E_Polynomial[i].power<<"\t"<<E_Polynomial[i].prob<<endl;
        for (j = 0; j < T_Polynomial.size(); j++) {
            temp_Polynomial.prob = E_Polynomial[i].prob * T_Polynomial[j].prob;
            if (temp_Polynomial.prob > MIN_PROB) {
                temp_Polynomial.power =
                    E_Polynomial[i].power + T_Polynomial[j].power;

                if (flag_size > FINE_GRID) {
                    index = unsigned(fabs(temp_Polynomial.power - min_mass) /
                                         delta_mass +
                                     0.5);
                    CGID_Polynomial[index].prob =
                        CGID_Polynomial[index].prob + temp_Polynomial.prob;
                    CGID_Polynomial[index].power =
                        CGID_Polynomial[index].power +
                        temp_Polynomial.prob * temp_Polynomial.power;
                } else {
                    CGID_Polynomial[index].prob =
                        CGID_Polynomial[index].prob + temp_Polynomial.prob;
                    CGID_Polynomial[index].power =
                        CGID_Polynomial[index].power +
                        temp_Polynomial.prob * temp_Polynomial.power;
                    index = index + 1;
                }
            }
        }
    }

    index = T_Polynomial.size();
    j = 0;
    for (i = 0; i < CGID_Polynomial.size(); i++) {
        if (CGID_Polynomial[i].prob != 0) {
            if (j < index) {
                T_Polynomial[j].prob = CGID_Polynomial[i].prob;
                T_Polynomial[j].power =
                    CGID_Polynomial[i].power / CGID_Polynomial[i].prob;
                j = j + 1;
            } else {
                temp_Polynomial.prob = CGID_Polynomial[i].prob;
                temp_Polynomial.power =
                    CGID_Polynomial[i].power / CGID_Polynomial[i].prob;
                T_Polynomial.push_back(temp_Polynomial);
            }
        }

        CGID_Polynomial[i].prob = 0;
        CGID_Polynomial[i].power = 0;
    }
    if (j < index) {
        T_Polynomial.erase(T_Polynomial.begin() + j, T_Polynomial.end());
    }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*Merging final coarse polynomial in Coarse Grid*/
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
vector<struct Polynomial> MIDAs::Merge_Coarse_Polynomial(
    vector<struct Polynomial> T_Polynomial) {
    int i, j, k;
    vector<struct Polynomial> Fine_Polynomial, F_Polynomial;
    struct Polynomial T;

    sort(T_Polynomial.begin(), T_Polynomial.end());

    for (k = 1; k <= 5; k++) {
        for (i = 0; i < int(T_Polynomial.size()) - 1; i++) {
            T.power = T_Polynomial[i].power * T_Polynomial[i].prob;
            T.prob = T_Polynomial[i].prob;
            if (T_Polynomial[i].power != 0) {
                for (j = i + 1; j < int(T_Polynomial.size()); j++) {
                    if (T_Polynomial[j].power != 0) {
                        if (fabs((T_Polynomial[i].power * MW_RESOLUTION -
                                  T_Polynomial[j].power * MW_RESOLUTION)) <=
                            k * MERGE_COARSE_RESOLUTION / 5) {
                            T.power = T.power + T_Polynomial[j].power *
                                                    T_Polynomial[j].prob;
                            T.prob = T.prob + T_Polynomial[j].prob;
                            T_Polynomial[i].power = T.power / T.prob;
                            T_Polynomial[i].prob = T.prob;
                            T_Polynomial[j].power = 0;
                            T_Polynomial[j].prob = 0;
                        } else {
                            break;
                        }
                    }
                }
                T_Polynomial[i].power = T.power / T.prob;
                T_Polynomial[i].prob = T.prob;
            }
        }
    }

    for (i = 0; i < int(T_Polynomial.size()); i++) {
        if (T_Polynomial[i].power != 0) {
            Fine_Polynomial.push_back(T_Polynomial[i]);
        }
    }

    return Fine_Polynomial;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*Merging final coarse polynomial in Coarse Grid*/
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
vector<struct Polynomial> MIDAs::Merge_Fine_Polynomial(
    vector<struct Polynomial> T_Polynomial) {
    int i, j, k;
    vector<struct Polynomial> Fine_Polynomial;
    struct Polynomial T;

    sort(T_Polynomial.begin(), T_Polynomial.end());
    for (k = 1; k <= 9; k++) {
        for (i = 0; i < int(T_Polynomial.size()) - 1; i++) {
            T.power = T_Polynomial[i].power * T_Polynomial[i].prob;
            T.prob = T_Polynomial[i].prob;
            if (T_Polynomial[i].power != 0) {
                for (j = i + 1; j < int(T_Polynomial.size()); j++) {
                    if (T_Polynomial[j].power != 0) {
                        if (k <= 8) {
                            if (fabs((T_Polynomial[i].power * MW_RESOLUTION -
                                      T_Polynomial[j].power * MW_RESOLUTION)) <=
                                (k * MERGE_FINE_RESOLUTION / 8)) {
                                T.power = T.power + T_Polynomial[j].power *
                                                        T_Polynomial[j].prob;
                                T.prob = T.prob + T_Polynomial[j].prob;
                                T_Polynomial[i].power = T.power / T.prob;
                                T_Polynomial[i].prob = T.prob;
                                T_Polynomial[j].power = 0;
                                T_Polynomial[j].prob = 0;
                            } else {
                                break;
                            }
                        } else {
                            if (fabs((T_Polynomial[i].power * MW_RESOLUTION -
                                      T_Polynomial[j].power * MW_RESOLUTION)) <=
                                (MERGE_FINE_RESOLUTION +
                                 MERGE_FINE_RESOLUTION / 100)) {
                                T.power = T.power + T_Polynomial[j].power *
                                                        T_Polynomial[j].prob;
                                T.prob = T.prob + T_Polynomial[j].prob;
                                T_Polynomial[i].power = T.power / T.prob;
                                T_Polynomial[i].prob = T.prob;
                                T_Polynomial[j].power = 0;
                                T_Polynomial[j].prob = 0;
                            } else {
                                break;
                            }
                        }
                    }
                }
                T_Polynomial[i].power = T.power / T.prob;
                T_Polynomial[i].prob = T.prob;
            }
        }
    }

    for (i = 0; i < int(T_Polynomial.size()); i++) {
        if (T_Polynomial[i].power != 0) {
            Fine_Polynomial.push_back(T_Polynomial[i]);
        }
    }

    return Fine_Polynomial;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*Computing FT-CGID resolution=1.0 dalton.
This algorithm is homologous to:
Femandez Cossio Diaz. Computation Of Isotopic peak Center-Mass Distriubtion by
Fourier Transform. Analytic Chemistry 2012, 84(16):7052-7056.*/
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void MIDAs::FT_JFC_ID(
    vector<vector<struct Composition> > &element_composition) {
    int i, k, j, N;
    double prob_min = 0;
    double mw_monoisotopic, np, sum_p;
    double mw, average_mass, mass_range, freq, delta;
    double resolution;
    double angle, radius, x, y, phi;
    double sx, sy;
    double dx, dy, nx, ny, angle_d, angle_n;
    double one_pi = acos(-1.0);
    double two_pi = 2 * one_pi;
    struct Isotopic_Distribution temp;
    vector<struct Isotopic_Distribution> B, A;
    cerr.precision(16);

    resolution = 1.0;
    // Computing mass range  2^n
    mass_range = ft_variance_molecular_mass(element_composition);
    mass_range = 16 * sqrt(1 + mass_range);
    mass_range = int(log(mass_range) / log(2.0) + 1.0);
    mass_range = int(pow(2.0, mass_range));
    N = int((mass_range / resolution));
    N = int(log(N) / log(2.0));
    N = int(pow(2.0, N));
    MASS_FREQUENCY_DOMAIN = (double *)malloc((2 * N + 1) * sizeof(double));
    MASS_FREQUENCY_DOMAIN_DP = (double *)malloc((2 * N + 1) * sizeof(double));
    delta = 1.0 / N;
    resolution = mass_range / N;

    average_mass = ft_average_molecular_mass(element_composition);
    mw_monoisotopic = monoisotopic_molecular_mass(element_composition);
    average_mass = int(average_mass - mw_monoisotopic + 0.5);
    cerr << "ID:---------------------------------------------------------------"
            "-------------------------"
         << endl;
    cerr << "ID:"
         << "monoisotopic mass = " << mw_monoisotopic << endl;
    cerr << "ID:"
         << "Average mass = " << ft_average_molecular_mass(element_composition)
         << endl;
    cerr << "ID:"
         << "FT resolution = " << resolution << endl;
    cerr << "ID:Number of grid points used = " << N << endl;
    cerr << "ID:Heterodyne mass = " << average_mass << endl;

    for (i = 1; i <= N / 2; i++) {
        freq = (i - 1) * delta;
        phi = angle = sx = sy = 0.0;
        radius = 1.0;
        for (k = 0; k < int(element_composition.size()); k++) {
            if (int(element_composition[k].size()) > 0) {
                x = y = 0;
                dx = dy = nx = ny = 0;
                for (j = 0; j < int(element_composition[k].size()); j++) {
                    // phase
                    mw = element_composition[k][j].nucleon / resolution;
                    phi = two_pi * mw * freq;

                    // probability convolution
                    x = x + element_composition[k][j].prob * cos(phi);
                    y = y + element_composition[k][j].prob * sin(phi);

                    // nucleon number
                    dx = dx + element_composition[k][j].prob * cos(phi);
                    dy = dy + element_composition[k][j].prob * sin(phi);

                    // average weight
                    nx = nx + element_composition[k][j].ave_mw * cos(phi);
                    ny = ny + element_composition[k][j].ave_mw * sin(phi);
                }

                radius = radius * pow(sqrt(x * x + y * y),
                                      element_composition[k][0].atoms);
                angle = angle + element_composition[k][0].atoms *
                                    atan2(y, x);  //[-pi,pi];

                angle_n = atan2(ny, nx);
                angle_d = atan2(dy, dx);
                sx = sx + (element_composition[k][0].atoms *
                           sqrt(nx * nx + ny * ny) / sqrt(dx * dx + dy * dy)) *
                              cos(angle_n - angle_d);
                sy = sy + (element_composition[k][0].atoms *
                           sqrt(nx * nx + ny * ny) / sqrt(dx * dx + dy * dy)) *
                              sin(angle_n - angle_d);
            }
        }
        // compute b(v)-term in JFC
        angle_n = atan2(sy, sx);
        x = sqrt(sx * sx + sy * sy) * radius *
            cos(angle + angle_n - two_pi * average_mass * freq);
        y = sqrt(sx * sx + sy * sy) * radius *
            sin(angle + angle_n - two_pi * average_mass * freq);

        MASS_FREQUENCY_DOMAIN[2 * i - 1] = x;
        MASS_FREQUENCY_DOMAIN[2 * i] = y;

        // compute probability a(v)-term in JFC
        x = radius * cos(angle - two_pi * average_mass * freq);
        y = radius * sin(angle - two_pi * average_mass * freq);

        MASS_FREQUENCY_DOMAIN_DP[2 * i - 1] = x;
        MASS_FREQUENCY_DOMAIN_DP[2 * i] = y;
    }

    for (i = N / 2 + 1; i <= N; i++) {
        freq = (i - N - 1) * delta;
        phi = angle = sx = sy = 0.0;
        radius = 1.0;
        for (k = 0; k < int(element_composition.size()); k++) {
            if (int(element_composition[k].size()) > 0) {
                x = y = 0;
                dx = dy = nx = ny = 0;
                for (j = 0; j < int(element_composition[k].size()); j++) {
                    // phase
                    mw = element_composition[k][j].nucleon / resolution;
                    phi = two_pi * mw * freq;

                    // probability convolution
                    x = x + element_composition[k][j].prob * cos(phi);
                    y = y + element_composition[k][j].prob * sin(phi);

                    // nucleon number
                    dx = dx + element_composition[k][j].prob * cos(phi);
                    dy = dy + element_composition[k][j].prob * sin(phi);

                    // average weight
                    nx = nx + element_composition[k][j].ave_mw * cos(phi);
                    ny = ny + element_composition[k][j].ave_mw * sin(phi);
                }

                radius = radius * pow(sqrt(x * x + y * y),
                                      element_composition[k][0].atoms);
                angle = angle + element_composition[k][0].atoms *
                                    atan2(y, x);  //[-pi,pi];

                angle_d = atan2(dy, dx);
                angle_n = atan2(ny, nx);
                sx = sx + (element_composition[k][0].atoms *
                           sqrt(nx * nx + ny * ny) / sqrt(dx * dx + dy * dy)) *
                              cos(angle_n - angle_d);
                sy = sy + (element_composition[k][0].atoms *
                           sqrt(nx * nx + ny * ny) / sqrt(dx * dx + dy * dy)) *
                              sin(angle_n - angle_d);
            }
        }

        // compute  b(v)-term in JFC
        angle_n = atan2(sy, sx);
        x = sqrt(sx * sx + sy * sy) * radius *
            cos(angle + angle_n - two_pi * average_mass * freq);
        y = sqrt(sx * sx + sy * sy) * radius *
            sin(angle + angle_n - two_pi * average_mass * freq);
        MASS_FREQUENCY_DOMAIN[2 * i - 1] = x;
        MASS_FREQUENCY_DOMAIN[2 * i] = y;

        // compute probability a(v)-term in JFC
        x = radius * cos(angle - two_pi * average_mass * freq);
        y = radius * sin(angle - two_pi * average_mass * freq);

        MASS_FREQUENCY_DOMAIN_DP[2 * i - 1] = x;
        MASS_FREQUENCY_DOMAIN_DP[2 * i] = y;
    }

    cerr << "ID:---------------------------------------------------------------"
            "-------------------------"
         << endl;

    // Performing inverse fourier transform. It normalizes the probabilities 1/N
    cerr << "ID:Start IDFFT" << endl;
    four1(MASS_FREQUENCY_DOMAIN, N, -1);
    four1(MASS_FREQUENCY_DOMAIN_DP, N, -1);
    cerr << "ID:Finished IDFFT" << endl;

    prob_min = 0;
    // probability used to filter oscillating terms
    for (i = 1; i <= N; i++) {
        y = MASS_FREQUENCY_DOMAIN[2 * i - 1] / N;
        if (y < prob_min) {
            prob_min = y;
        }
    }
    prob_min = 1e-12;  // probability used to filter oscillating terms
    cerr << "Mim prob=" << prob_min << endl;

    // selecting terms above min probability
    np = 0;
    for (i = N / 2 + 1; i <= N; i++) {
        temp.mw = (i - N - 1) * resolution + average_mass * resolution;
        temp.prob = MASS_FREQUENCY_DOMAIN[2 * i - 1] / N;
        // cout<<temp.mw<<"\t"<<temp.prob<<endl;

        if (temp.prob > prob_min &&
            MASS_FREQUENCY_DOMAIN_DP[2 * i - 1] / N > prob_min) {
            B.push_back(temp);
            temp.prob = MASS_FREQUENCY_DOMAIN_DP[2 * i - 1] / N;
            A.push_back(temp);
            np = np + temp.prob;
        }
    }
    for (i = 1; i <= N / 2 - 1; i++) {
        temp.mw = (i - 1) * resolution + average_mass * resolution;
        temp.prob = MASS_FREQUENCY_DOMAIN[2 * i - 1] / N;

        if (temp.prob > prob_min &&
            MASS_FREQUENCY_DOMAIN_DP[2 * i - 1] / N > prob_min) {
            B.push_back(temp);
            temp.prob = MASS_FREQUENCY_DOMAIN_DP[2 * i - 1] / N;
            A.push_back(temp);
            np = np + temp.prob;
        }
    }

    // normalizing probability
    x = sum_p = 0;
    for (i = 0; i < int(B.size()); i++) {
        // A[i].prob=A[i].prob/np;
        A[i].prob = A[i].prob;
        sum_p = sum_p + A[i].prob;

        B[i].mw = B[i].prob / A[i].prob;  // computes average molecular weight
        B[i].prob = A[i].prob;            // computes probability
        x = x + B[i].mw * B[i].prob;
        // cout<<B[i].mw<<"\t"<<B[i].prob<<endl;
        COARSE_GRAINED_ISOTOPIC_DISTRIBUTION.push_back(B[i]);
    }

    cerr << "ID:Sum of probabilities CGID = " << sum_p << endl;
    cerr << "ID:Computed average mass = " << x << "\t" << x / sum_p << endl;

    cerr << "ID:---------------------------------------------------------------"
            "-------------------------"
         << endl;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*Computing FT-CGID with low resolution.
This algorithm is homologous to:
Alan L. Rockwood. Rapid Calculation of Isotope Distriubtion.
Analytical Chemistry, 1995, 67, 2699-2704. (Mercury)*/
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void MIDAs::FT_Coarse_Grained_ID(
    vector<vector<struct Composition> > &element_composition,
    vector<struct Isotopic_Distribution> &COARSE_GRAINED_ISOTOPIC_DISTRIBUTION,
    double resolution) {
    int i, k, j, flag, N;
    double prob_min = 0;
    double mw_monoisotopic, np, sum_p;
    double ave1, ave2, sigma1, sigma2, ratio;
    double mw, average_mass, mass_range, freq, delta;
    double rho_resolution, used_resolution;
    double angle, radius, x, y, phi;
    double one_pi = acos(-1.0);
    double two_pi = 2 * one_pi;
    struct Isotopic_Distribution temp;
    vector<struct Isotopic_Distribution> T_ID;
    cerr.precision(12);

    // Computing mass range  2^n
    used_resolution = resolution;
    rho_resolution = 1.0;
    delta = 1.0;
    flag = N = k = 0;
    while (flag == 0) {
        resolution = used_resolution / pow(2.0, k);
        mass_range = ft_variance_molecular_mass(element_composition);
        mass_range = int(15.0 * sqrt(1 + mass_range) + 1);
        mass_range = int(log(mass_range) / log(2.0) + 1.0);
        mass_range = int(pow(2.0, mass_range));
        N = int(mass_range / resolution);
        N = int((mass_range / resolution));
        N = int(log(N) / log(2.0) + 1);
        N = int(pow(2.0, N));
        delta = 1.0 / N;
        resolution = mass_range / N;
        rho_resolution = resolution;
        k = k + 1;
        if (rho_resolution < used_resolution) {
            MASS_FREQUENCY_DOMAIN =
                (double *)malloc((2 * N + 1) * sizeof(double));
            flag = 1;
        }
    }

    average_mass =
        int(ft_average_molecular_mass(element_composition, resolution) + 0.5);
    average_mass = average_mass / resolution;
    mw_monoisotopic = monoisotopic_molecular_mass(element_composition);
    cerr << "ID:---------------------------------------------------------------"
            "-------------------------"
         << endl;
    cerr << "ID:"
         << "Lightest mass = " << mw_monoisotopic << endl;
    cerr << "ID:"
         << "Average mass = "
         << ft_average_molecular_mass(element_composition, resolution) << endl;
    cerr << "ID:"
         << "FT resolution = " << resolution << endl;
    cerr << "ID:Number of grid points used = " << N << endl;
    cerr << "ID:Heterodyne mass = " << average_mass << endl;

    for (i = 1; i <= N / 2; i++) {
        freq = (i - 1) * delta;

        phi = angle = 0.0;
        radius = 1.0;
        for (k = 0; k < int(element_composition.size()); k++) {
            if (int(element_composition[k].size()) > 0) {
                x = y = 0;
                for (j = 0; j < int(element_composition[k].size()); j++) {
                    mw = int(element_composition[k][j].mw / resolution + 0.5);
                    phi = two_pi * mw * freq;
                    x = x + element_composition[k][j].prob * cos(phi);
                    y = y + element_composition[k][j].prob * sin(phi);
                }

                radius = radius * pow(sqrt(x * x + y * y),
                                      element_composition[k][0].atoms);
                angle = angle + element_composition[k][0].atoms *
                                    atan2(y, x);  //[-pi,pi];
            }
        }
        x = radius * cos(angle - two_pi * average_mass * freq);
        y = radius * sin(angle - two_pi * average_mass * freq);

        MASS_FREQUENCY_DOMAIN[2 * i - 1] = x;
        MASS_FREQUENCY_DOMAIN[2 * i] = y;
    }

    for (i = N / 2 + 1; i <= N; i++) {
        freq = (i - N - 1) * delta;

        phi = angle = 0.0;
        radius = 1.0;
        for (k = 0; k < int(element_composition.size()); k++) {
            if (int(element_composition[k].size()) > 0) {
                x = y = 0;
                for (j = 0; j < int(element_composition[k].size()); j++) {
                    mw = int(element_composition[k][j].mw / resolution + 0.5);
                    phi = two_pi * mw * freq;
                    x = x + element_composition[k][j].prob * cos(phi);
                    y = y + element_composition[k][j].prob * sin(phi);
                }

                radius = radius * pow(sqrt(x * x + y * y),
                                      element_composition[k][0].atoms);
                angle = angle + element_composition[k][0].atoms *
                                    atan2(y, x);  //[-pi,pi];
            }
        }
        x = radius * cos(angle - two_pi * average_mass * freq);
        y = radius * sin(angle - two_pi * average_mass * freq);

        MASS_FREQUENCY_DOMAIN[2 * i - 1] = x;
        MASS_FREQUENCY_DOMAIN[2 * i] = y;
    }
    cerr << "ID:---------------------------------------------------------------"
            "-------------------------"
         << endl;

    // Performing inverse fourier transform. It normalizes the probabilities 1/N
    cerr << "ID:Start IDFFT" << endl;
    four1(MASS_FREQUENCY_DOMAIN, N, -1);
    cerr << "ID:Finished IDFFT" << endl;

    // Correcting mass and removing oscillating terms
    ave1 = ft_average_molecular_mass(element_composition);
    ave2 = ft_average_molecular_mass(element_composition, resolution);
    // ave3 = int(ave2 + 0.5);
    sigma1 = sqrt(ft_variance_molecular_mass(element_composition));
    sigma2 = sqrt(ft_variance_molecular_mass(element_composition, resolution));
    ratio = sigma1 / sigma2;

    // probability used to filter oscillating terms
    prob_min = 0;
    for (i = 1; i <= N; i++) {
        y = MASS_FREQUENCY_DOMAIN[2 * i - 1] / N;
        if (y < prob_min) {
            prob_min = y;
        }
    }
    prob_min = -2 * prob_min;  // probability used to filter oscillating terms
    cerr << "ID:CGID ratio =" << ratio << endl;
    cerr << "ID:CGID leakage probability =" << prob_min << endl;

    // selecting terms above min probability
    for (i = N / 2 + 1; i <= N; i++) {
        temp.mw = ratio * ((i - N - 1) * resolution +
                           average_mass * resolution - ave2) +
                  ave1;
        temp.prob = MASS_FREQUENCY_DOMAIN[2 * i - 1] / N;

        if (temp.prob > prob_min) {
            T_ID.push_back(temp);
        }
    }
    for (i = 1; i <= N / 2 - 1; i++) {
        temp.mw =
            ratio * ((i - 1) * resolution + average_mass * resolution - ave2) +
            ave1;
        temp.prob = MASS_FREQUENCY_DOMAIN[2 * i - 1] / N;

        if (temp.prob > prob_min) {
            T_ID.push_back(temp);
        }
    }

    // Merging terms terms within requested resolution
    for (k = 1; k <= 5; k++) {
        for (i = 0; i < int(T_ID.size()) - 1; i++) {
            if (T_ID[i].mw > (mw_monoisotopic - 0.01)) {
                temp.mw = T_ID[i].mw * T_ID[i].prob;
                temp.prob = T_ID[i].prob;
                if (T_ID[i].mw != 0) {
                    for (j = i + 1; j < int(T_ID.size()); j++) {
                        if (T_ID[j].mw != 0) {
                            if (fabs(T_ID[i].mw - T_ID[j].mw) <=
                                k * MERGE_COARSE_RESOLUTION / 5.0) {
                                temp.mw = temp.mw + T_ID[j].mw * T_ID[j].prob;
                                temp.prob = temp.prob + T_ID[j].prob;
                                T_ID[i].mw = temp.mw / temp.prob;
                                T_ID[i].prob = temp.prob;
                                T_ID[j].mw = 0;
                                T_ID[j].prob = 0;
                            } else {
                                break;
                            }
                        }
                    }
                    T_ID[i].mw = temp.mw / temp.prob;
                    T_ID[i].prob = temp.prob;
                }
            } else {
                T_ID[i].mw = 0;
                T_ID[i].prob = 0;
            }
        }
    }

    // Normalizing CGID probability
    sum_p = np = 0;
    for (i = 0; i < int(T_ID.size()); i++) {
        if (T_ID[i].mw > 0) {
            np = np + T_ID[i].prob;
        }
    }
    cerr << "ID:Sum of probabilities CGID = " << np << endl;

    for (i = 0; i < int(T_ID.size()); i++) {
        if (T_ID[i].mw > 0) {
            T_ID[i].prob = T_ID[i].prob / np;
            sum_p = sum_p + T_ID[i].prob;
            COARSE_GRAINED_ISOTOPIC_DISTRIBUTION.push_back(T_ID[i]);
        }
    }
    cerr << "ID:Sum of probabilities CGID = " << sum_p << endl;
    cerr << "ID:---------------------------------------------------------------"
            "-------------------------"
         << endl;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*Computing FT-FGID with adjustable resolution.
This algorithm is homologous to:
Alan L. Rockwood. Rapid Calculation of Isotope Distriubtion.
Analytical Chemistry, 1995, 67, 2699-2704. (Mercury)*/
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void MIDAs::FT_Fine_Grained_ID(
    vector<vector<struct Composition> > &element_composition,
    vector<struct Polynomial> &FINE_ISOTOPIC_DISTRIBUTION, double resolution) {
    int i, k, j, flag;
    long int N;
    double prob_min = 0;
    double mw_monoisotopic, np, sum_p;
    double ave1, ave2, sigma1, sigma2, ratio;
    double mw, average_mass, mass_range, freq, delta;
    double average_prob = 0;
    double rho_resolution, used_resolution;
    double angle, radius, x, y, phi;
    double one_pi = acos(-1);
    double two_pi = 2 * one_pi;
    struct Polynomial temp;
    vector<struct Polynomial> T_ID;
    cerr.precision(12);

    // Computing mass range  2^n
    used_resolution = resolution;
    rho_resolution = 1.0;
    delta = 1.0;
    flag = N = k = 0;
    while (flag == 0) {
        resolution = used_resolution / pow(2.0, k);
        mass_range = ft_variance_molecular_mass(element_composition);
        mass_range = int(15.0 * sqrt(1 + mass_range) + 1);
        mass_range = int(log(mass_range) / log(2.0) + 1.0);
        mass_range = int(pow(2.0, mass_range));
        N = int(mass_range / resolution);
        N = int(log(N) / log(2.0) + 1.0);
        N = int(pow(2.0, N));
        delta = 1.0 / N;
        resolution = mass_range / N;
        rho_resolution = resolution;
        k = k + 1;
        if (rho_resolution < used_resolution) {
            MASS_FREQUENCY_DOMAIN =
                (double *)malloc((2 * N + 1) * sizeof(double));
            flag = 1;
        }
    }

    average_mass =
        int(ft_average_molecular_mass(element_composition, resolution) + 0.5);
    average_mass = average_mass / resolution;
    mw_monoisotopic = monoisotopic_molecular_mass(element_composition);
    cerr << "ID:---------------------------------------------------------------"
            "-------------------------"
         << endl;
    cerr << "ID:"
         << "Lightest mass = " << mw_monoisotopic << endl;
    cerr << "ID:"
         << "Average mass = "
         << ft_average_molecular_mass(element_composition, resolution) << endl;
    cerr << "ID:"
         << "FT resolution = " << resolution << endl;
    cerr << "ID:Number of grid points used = " << N << endl;
    cerr << "ID:Heterodyne mass = " << average_mass << endl;

    for (i = 1; i <= N / 2; i++) {
        freq = (i - 1) * delta;

        phi = angle = 0.0;
        radius = 1.0;
        for (k = 0; k < int(element_composition.size()); k++) {
            if (int(element_composition[k].size()) > 0) {
                x = y = 0;
                for (j = 0; j < int(element_composition[k].size()); j++) {
                    mw = int(element_composition[k][j].mw / resolution + 0.5);
                    phi = two_pi * mw * freq;
                    x = x + element_composition[k][j].prob * cos(phi);
                    y = y + element_composition[k][j].prob * sin(phi);
                }

                radius = radius * pow(sqrt(x * x + y * y),
                                      element_composition[k][0].atoms);
                angle = angle + element_composition[k][0].atoms *
                                    atan2(y, x);  //[-pi,pi];
            }
        }
        x = radius * cos(angle - two_pi * average_mass * freq);
        y = radius * sin(angle - two_pi * average_mass * freq);

        MASS_FREQUENCY_DOMAIN[2 * i - 1] = x;
        MASS_FREQUENCY_DOMAIN[2 * i] = y;
    }

    for (i = N / 2 + 1; i <= N; i++) {
        freq = (i - N - 1) * delta;

        phi = angle = 0.0;
        radius = 1.0;
        for (k = 0; k < int(element_composition.size()); k++) {
            if (int(element_composition[k].size()) > 0) {
                x = y = 0;
                for (j = 0; j < int(element_composition[k].size()); j++) {
                    mw = int(element_composition[k][j].mw / resolution + 0.5);
                    phi = two_pi * mw * freq;
                    x = x + element_composition[k][j].prob * cos(phi);
                    y = y + element_composition[k][j].prob * sin(phi);
                }

                radius = radius * pow(sqrt(x * x + y * y),
                                      element_composition[k][0].atoms);
                angle = angle + element_composition[k][0].atoms *
                                    atan2(y, x);  //[-pi,pi];
            }
        }
        x = radius * cos(angle - two_pi * average_mass * freq);
        y = radius * sin(angle - two_pi * average_mass * freq);

        MASS_FREQUENCY_DOMAIN[2 * i - 1] = x;
        MASS_FREQUENCY_DOMAIN[2 * i] = y;
    }
    cerr << "ID:---------------------------------------------------------------"
            "-------------------------"
         << endl;
    cerr << "ID:Start IDFFT" << endl;
    // Performing inverse fourier transform. It normalizes the probabilities 1/N
    four1(MASS_FREQUENCY_DOMAIN, N, -1);
    cerr << "ID:Finished IDFFT" << endl;

    // Correcting mass and removing oscillating terms
    ave1 = ft_average_molecular_mass(element_composition);
    ave2 = ft_average_molecular_mass(element_composition, resolution);
    // ave3 = int(ave2 + 0.5);
    sigma1 = sqrt(ft_variance_molecular_mass(element_composition));
    sigma2 = sqrt(ft_variance_molecular_mass(element_composition, resolution));
    ratio = sigma1 / sigma2;

    // probability used to filter oscillating terms
    prob_min = 0;
    for (i = 1; i <= N; i++) {
        y = MASS_FREQUENCY_DOMAIN[2 * i - 1] / N;
        if (y < prob_min) {
            prob_min = y;
        }
    }
    prob_min = -2.0 * prob_min;

    k = 0;
    for (i = 1; i <= N; i++) {
        y = MASS_FREQUENCY_DOMAIN[2 * i - 1] / N;
        if (y > prob_min) {
            average_prob = average_prob + y;
            k = k + 1;
        }
    }
    if (k > 0) {
        average_prob = average_prob / k;
    }

    cerr << "ID:FGID ratio =" << ratio << endl;
    cerr << "ID:FGID leakage probability =" << prob_min << endl;
    cerr << "ID:FGID average_prob =" << average_prob << endl;

    // Filtering probability terms
    for (i = N / 2 + 1; i <= N; i++) {
        temp.power = ratio * ((i - N - 1) * resolution +
                              average_mass * resolution - ave2) +
                     ave1;
        temp.prob = MASS_FREQUENCY_DOMAIN[2 * i - 1] / N;

        if (temp.prob > prob_min) {
            T_ID.push_back(temp);
        }
    }
    for (i = 1; i <= N / 2 - 1; i++) {
        temp.power =
            ratio * ((i - 1) * resolution + average_mass * resolution - ave2) +
            ave1;
        temp.prob = MASS_FREQUENCY_DOMAIN[2 * i - 1] / N;

        if (temp.prob > prob_min) {
            T_ID.push_back(temp);
        }
    }

    // Merging terms within resolution
    for (k = 1; k <= 9; k++) {
        for (i = 0; i < int(T_ID.size()) - 1; i++) {
            if (T_ID[i].power >= mw_monoisotopic) {
                temp.power = T_ID[i].power * T_ID[i].prob;
                temp.prob = T_ID[i].prob;
                if (T_ID[i].power != 0) {
                    for (j = i + 1; j < int(T_ID.size()); j++) {
                        if (T_ID[j].power != 0) {
                            if (k <= 8) {
                                if (fabs(T_ID[i].power - T_ID[j].power) <=
                                    (k * MERGE_FINE_RESOLUTION / 8)) {
                                    temp.power = temp.power +
                                                 T_ID[j].power * T_ID[j].prob;
                                    temp.prob = temp.prob + T_ID[j].prob;
                                    T_ID[i].power = temp.power / temp.prob;
                                    T_ID[i].prob = temp.prob;
                                    T_ID[j].power = 0;
                                    T_ID[j].prob = 0;
                                } else {
                                    break;
                                }
                            } else {
                                if (fabs(T_ID[i].power - T_ID[j].power) <=
                                    (MERGE_FINE_RESOLUTION +
                                     MERGE_FINE_RESOLUTION / 100)) {
                                    temp.power = temp.power +
                                                 T_ID[j].power * T_ID[j].prob;
                                    temp.prob = temp.prob + T_ID[j].prob;
                                    T_ID[i].power = temp.power / temp.prob;
                                    T_ID[i].prob = temp.prob;
                                    T_ID[j].power = 0;
                                    T_ID[j].prob = 0;
                                } else {
                                    break;
                                }
                            }
                        }
                    }
                    T_ID[i].power = temp.power / temp.prob;
                    T_ID[i].prob = temp.prob;
                } else {
                    T_ID[i].power = 0;
                    T_ID[i].prob = 0;
                }
            }
        }
    }

    // Normalizing CGID probability
    sum_p = np = 0;
    for (i = 0; i < int(T_ID.size()); i++) {
        if (T_ID[i].power > 0) {
            np = np + T_ID[i].prob;
        }
    }
    cerr << "ID:Sum of probabilities FGID = " << np << endl;

    for (i = 0; i < int(T_ID.size()); i++) {
        if (T_ID[i].power > 0) {
            T_ID[i].prob = T_ID[i].prob / np;
            sum_p = sum_p + T_ID[i].prob;
            FINE_ISOTOPIC_DISTRIBUTION.push_back(T_ID[i]);
        }
    }
    cerr << "ID:Sum of probabilities FGID = " << sum_p << endl;

    cerr << "ID:---------------------------------------------------------------"
            "-------------------------"
         << endl;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*The function below computes the fast fourier transform. 1==FFT,-1=IFFT.*/
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void MIDAs::four1(double *Data, int nn, int isign) {
    unsigned long i, j, m, n, mmax, istep;
    double wr, wpr, wpi, wi, theta;
    double wtemp, tempr, tempi;
    double one_pi = acos(-1);
    double two_pi = 2 * one_pi;

    /* Perform bit reversal of Data[] */
    n = nn << 1;
    j = 1;
    for (i = 1; i < n; i += 2) {
        if (j > i) {
            wtemp = Data[i];
            Data[i] = Data[j];
            Data[j] = wtemp;
            wtemp = Data[i + 1];
            Data[i + 1] = Data[j + 1];
            Data[j + 1] = wtemp;
        }
        m = n >> 1;
        while (m >= 2 && j > m) {
            j -= m;
            m >>= 1;
        }
        j += m;
    }

    /* Perform Danielson-Lanczos section of FFT */
    n = nn << 1;
    mmax = 2;
    while (n > mmax) /* Loop executed log(2)nn times */
    {
        istep = mmax << 1;
        theta = isign * (two_pi / mmax);
        wtemp = sin(0.5 * theta);
        wpr = -2.0 * wtemp * wtemp;
        wpi = sin(theta);
        wr = 1.0;
        wi = 0.0;
        for (m = 1; m < mmax; m += 2) {
            for (i = m; i <= n; i += istep) {
                j = i + mmax;

                tempr = wr * Data[j] - wi * Data[j + 1];
                tempi = wr * Data[j + 1] + wi * Data[j];
                Data[j] = Data[i] - tempr;
                Data[j + 1] = Data[i + 1] - tempi;
                Data[i] += tempr;
                Data[i + 1] += tempi;
            }
            wr = (wtemp = wr) * wpr - wi * wpi + wr;
            wi = wi * wpr + wtemp * wpi + wi;
        }
        mmax = istep;
    }
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*Computes element average molecular weight for FT*/
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double MIDAs::ft_average_molecular_mass(
    vector<vector<struct Composition> > &element_composition,
    double resolution) {
    int k, i;
    double mw, ave_mw;

    mw = ave_mw = 0.0;

    for (k = 0; k < int(element_composition.size()); k++) {
        for (i = 0; i < int(element_composition[k].size()); i++) {
            if (resolution == 0) {
                ave_mw = ave_mw + element_composition[k][i].atoms *
                                      element_composition[k][i].mw *
                                      element_composition[k][i].prob;
            } else {
                mw = int(element_composition[k][i].mw / resolution + 0.5);
                mw = mw * resolution;

                ave_mw = ave_mw + element_composition[k][i].atoms * mw *
                                      element_composition[k][i].prob;
            }
        }
    }

    return ave_mw;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*Computes element variances. Returns theoretical variance for FT*/
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double MIDAs::ft_variance_molecular_mass(
    vector<vector<struct Composition> > &element_composition,
    double resolution) {
    int k, i;
    double sigma, var_mw, ave_mw, mw;

    mw = sigma = var_mw = ave_mw = 0;
    for (k = 0; k < int(element_composition.size()); k++) {
        var_mw = ave_mw = 0;
        if (int(element_composition[k].size()) > 0) {
            for (i = 0; i < int(element_composition[k].size()); i++) {
                if (resolution == 0) {
                    ave_mw = ave_mw + element_composition[k][i].mw *
                                          element_composition[k][i].prob;
                } else {
                    mw = int(element_composition[k][i].mw / resolution + 0.5);
                    mw = mw * resolution;
                    ave_mw = ave_mw + mw * element_composition[k][i].prob;
                }
            }

            for (i = 0; i < int(element_composition[k].size()); i++) {
                if (resolution == 0) {
                    var_mw =
                        var_mw +
                        element_composition[k][i].prob *
                            pow(element_composition[k][i].mw - ave_mw, 2.0);
                } else {
                    mw = int(element_composition[k][i].mw / resolution + 0.5);
                    mw = mw * resolution;
                    var_mw = var_mw + element_composition[k][i].prob *
                                          pow(mw - ave_mw, 2.0);
                }
            }

            var_mw = element_composition[k][0].atoms * var_mw;
            sigma = sigma + var_mw;
        }
    }
    return sigma;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*Sort element composition by increasing mw and normalize probability*/
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void MIDAs::sort_element_composition(
    vector<vector<struct Composition> > &Element_Composition) {
    int i, j, k;
    double sum_prob;
    struct Composition temp_Composition;

    for (k = 0; k < int(Element_Composition.size()); k++) {
        // Normalizing probability
        sum_prob = 0;
        for (i = 0; i < int(Element_Composition[k].size()); i++) {
            sum_prob = sum_prob + Element_Composition[k][i].prob;
        }
        for (i = 0; i < int(Element_Composition[k].size()); i++) {
            Element_Composition[k][i].prob =
                Element_Composition[k][i].prob / sum_prob;
        }

        sum_prob = 0;
        for (i = 0; i < int(Element_Composition[k].size()); i++) {
            Element_Composition[k][i].log_prob =
                log(Element_Composition[k][i].prob);
            Element_Composition[k][i].prob = Element_Composition[k][i].prob;
            Element_Composition[k][i].power =
                floor(Element_Composition[k][i].mw / MW_RESOLUTION + 0.5);
            sum_prob = sum_prob + Element_Composition[k][i].prob;
            cerr.precision(15);
        }

        if (Element_Composition[k].size() > 0) {
            if (fabs(sum_prob - 1.0) > 1e-15) {
                cerr.precision(16);
                cerr << "ID:Isotopic probabilities do not sum to 1. Sum = "
                     << sum_prob << endl;
                exit(1);
            }
        }
    }

    k = Element_Composition.size();
    for (k = 0; k < int(Element_Composition.size()); k++) {
        for (i = 0; i < int(Element_Composition[k].size()) - 1; i++) {
            for (j = i + 1; j < int(Element_Composition[k].size()); j++) {
                if (Element_Composition[k][i].mw >
                    Element_Composition[k][j].mw) {
                    temp_Composition = Element_Composition[k][i];
                    Element_Composition[k][i] = Element_Composition[k][j];
                    Element_Composition[k][j] = temp_Composition;
                }
            }
        }
    }
    // assigning nucleon number and average mass
    for (k = 0; k < int(Element_Composition.size()); k++) {
        for (i = 0; i < int(Element_Composition[k].size()); i++) {
            if (i == 0)  // lightest
            {
                Element_Composition[k][i].nucleon = 0;
                Element_Composition[k][i].ave_mw =
                    Element_Composition[k][i].mw *
                    Element_Composition[k][i].prob;
            } else {
                Element_Composition[k][i].nucleon =
                    int(Element_Composition[k][i].mw) -
                    int(Element_Composition[k][0].mw);
                Element_Composition[k][i].ave_mw =
                    Element_Composition[k][i].mw *
                    Element_Composition[k][i].prob;
            }
            cerr << "ID:" << Element_Composition[k][i].element << "\t"
                 << Element_Composition[k][i].prob << "\t"
                 << Element_Composition[k][i].mw << "\t"
                 << Element_Composition[k][i].atoms << "\t"
                 << Element_Composition[k][i].power << "\t"
                 << Element_Composition[k][i].ave_mw << "\t"
                 << Element_Composition[k][i].nucleon << "\t" << k << endl;
        }
    }

    cerr << "ID:Molecular element composition" << endl;
    cerr << "ID:element,probability, mass, number of atoms, adjusted mass by "
            "resolution, average mass, nucleon number"
         << endl;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*Adding charge state to molecule*/
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void MIDAs::add_charge_state(
    vector<vector<struct Composition> > &Element_Composition,
    int CHARGE_STATE) {
    int i, j, k;
    struct Composition temp_Composition;

    if (CHARGE_STATE > 0)  // Adding mass of hydrogen ion to be convoluted
    {
        Element_Composition.push_back(vector<struct Composition>());
        j = Element_Composition.size() - 1;
        for (i = 0; i < int(Element_Composition.size()); i++) {
            for (k = 0; k < int(Element_Composition[i].size()); k++) {
                if (strcmp(Element_Composition[i][k].element, "H") == 0) {
                    strcpy(temp_Composition.element, "H+");
                    temp_Composition.mw =
                        Element_Composition[i][k].mw - MASS_ELECTRON;
                    temp_Composition.prob = Element_Composition[i][k].prob;
                    temp_Composition.log_prob =
                        log(Element_Composition[i][k].prob);
                    temp_Composition.atoms = CHARGE_STATE;
                    Element_Composition[j].push_back(temp_Composition);
                }
            }
        }
    } else if (CHARGE_STATE < 0) {
        Element_Composition.push_back(vector<struct Composition>());
        j = Element_Composition.size() - 1;
        strcpy(temp_Composition.element, "e-");
        temp_Composition.mw = MASS_ELECTRON;
        temp_Composition.prob = 100.0;
        temp_Composition.log_prob = log(100.0);
        temp_Composition.atoms = fabs(CHARGE_STATE);
        Element_Composition[j].push_back(temp_Composition);
    }
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*Computes elements moments*/
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void MIDAs::Elements_Cumulants(vector<vector<double> > &cumulants, int n) {
    int i, j, l;
    double k;
    vector<vector<double> > moments;

    l = int(cumulants.size());
    for (i = 0; i < l; i++) {
        moments.push_back(vector<double>());
        for (j = 0; j < n; j++) {
            moments[i].push_back(0.0);
        }
    }

    // compute elements moments
    Elements_Raw_Moments(moments, n);

    // compute distribution cumulant
    for (j = 0; j < l; j++) {
        for (i = 1; i <= n; i++) {
            if (i == 1) {
                k = moments[j][0];

                cumulants[j][0] = k;
            } else if (i == 2) {
                k = moments[j][1] - pow(moments[j][0], 2.0);

                cumulants[j][1] = k;
            } else if (i == 3) {
                k = moments[j][2] - 3 * moments[j][1] * moments[j][0] +
                    2 * pow(moments[j][0], 3.0);

                cumulants[j][2] = k;
            } else if (i == 4) {
                k = moments[j][3] - 4 * moments[j][2] * moments[j][0] -
                    3 * pow(moments[j][1], 2.0) +
                    12 * moments[j][1] * pow(moments[j][0], 2.0) -
                    6 * pow(moments[j][0], 4.0);

                cumulants[j][3] = k;
            } else if (i == 5) {
                k = moments[j][4] - 5 * moments[j][3] * moments[j][0] -
                    10 * moments[j][2] * moments[j][1] +
                    20 * moments[j][2] * pow(moments[j][0], 2.0) +
                    30 * pow(moments[j][1], 2.0) * moments[j][0] -
                    60 * moments[j][1] * pow(moments[j][0], 3.0) +
                    24 * pow(moments[j][0], 5.0);

                cumulants[j][4] = k;
            } else if (i == 6) {
                k = moments[j][5] - 6 * moments[j][4] * moments[j][0] -
                    15 * moments[j][3] * moments[j][1] +
                    30 * moments[j][3] * pow(moments[j][0], 2.0) -
                    10 * pow(moments[j][2], 2.0) +
                    120 * moments[j][2] * moments[j][1] * moments[j][0] -
                    120 * moments[j][2] * pow(moments[j][0], 3.0) +
                    30 * pow(moments[j][1], 3.0) -
                    270 * pow(moments[j][1], 2.0) * pow(moments[j][0], 2.0) +
                    360 * moments[j][1] * pow(moments[j][0], 4.0) -
                    120 * pow(moments[j][0], 6.0);

                cumulants[j][5] = k;
            } else if (i == 7) {
                k = moments[j][6] - 7 * moments[j][5] * moments[j][0] -
                    21 * moments[j][4] * moments[j][1] +
                    42 * moments[j][4] * pow(moments[j][0], 2.0) -
                    35 * moments[j][3] * moments[j][2] +
                    210 * moments[j][3] * moments[j][1] * moments[j][0] -
                    210 * moments[j][3] * pow(moments[j][0], 3.0) +
                    140 * pow(moments[j][2], 2.0) * moments[j][0] +
                    210 * moments[j][2] * pow(moments[j][1], 2.0) -
                    1260 * moments[j][2] * moments[j][1] *
                        pow(moments[j][0], 2.0) +
                    840 * moments[j][2] * pow(moments[j][0], 4.0) -
                    630 * pow(moments[j][1], 3.0) * moments[j][0] +
                    2520 * pow(moments[j][1], 2.0) * pow(moments[j][0], 3.0) -
                    2520 * moments[j][1] * pow(moments[j][0], 5.0) +
                    720 * pow(moments[j][0], 7.0);

                cumulants[j][6] = k;
            } else if (i == 8) {
                k = moments[j][7] - 8 * moments[j][6] * moments[j][0] -
                    28 * moments[j][5] * moments[j][1] +
                    56 * moments[j][5] * pow(moments[j][0], 2.0) -
                    56 * moments[j][4] * moments[j][2] +
                    336 * moments[j][4] * moments[j][1] * moments[j][0] -
                    336 * moments[j][4] * pow(moments[j][0], 3.0) -
                    35 * pow(moments[j][3], 2.0) +
                    560 * moments[j][3] * moments[j][2] * moments[j][0] +
                    420 * moments[j][3] * pow(moments[j][1], 2.0) -
                    2520 * moments[j][3] * moments[j][1] *
                        pow(moments[j][0], 2.0) +
                    1680 * moments[j][3] * pow(moments[j][0], 4.0) +
                    560 * pow(moments[j][2], 2.0) * moments[j][1] -
                    1680 * pow(moments[j][2], 2.0) * pow(moments[j][0], 2.0) -
                    5040 * moments[j][2] * pow(moments[j][1], 2.0) *
                        moments[j][0] +
                    13440 * moments[j][2] * moments[j][1] *
                        pow(moments[j][0], 3.0) -
                    6720 * moments[j][2] * pow(moments[j][0], 5.0) -
                    630 * pow(moments[j][1], 4.0) +
                    10080 * pow(moments[j][1], 3.0) * pow(moments[j][0], 2.0) -
                    25200 * pow(moments[j][1], 2.0) * pow(moments[j][0], 4.0) +
                    20160 * moments[j][1] * pow(moments[j][0], 6.0) -
                    5040 * pow(moments[j][0], 8.0);

                cumulants[j][7] = k;
            }
        }
    }
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*Computes elements moments*/
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void MIDAs::Elements_Raw_Moments(vector<vector<double> > &elements_moments_mass,
                                 int n) {
    int i, j, k;
    double mw;
    for (k = 0; k < int(Element_Composition.size()); k++) {
        for (j = 0; j < n; j++) {
            mw = 0;
            for (i = 0; i < int(Element_Composition[k].size()); i++) {
                mw = mw + Element_Composition[k][i].prob *
                              pow(Element_Composition[k][i].mw, j + 1);
            }
            elements_moments_mass[k][j] = mw;
        }
    }
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*Return elements average mass*/
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void MIDAs::Elements_Average_Mass(vector<double> &elements_average_mass) {
    int k, i;
    double mw = 0.0;

    for (k = 0; k < int(Element_Composition.size()); k++) {
        mw = 0;
        for (i = 0; i < int(Element_Composition[k].size()); i++) {
            mw = mw +
                 Element_Composition[k][i].prob * Element_Composition[k][i].mw;
        }
        elements_average_mass.push_back(mw);
    }
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*Return elements average mass*/
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void MIDAs::Elements_Variance_Mass(vector<double> &elements_variance_mass) {
    int k, i;
    double v_mw = 0.0;
    double x, xx;
    for (k = 0; k < int(Element_Composition.size()); k++) {
        x = xx = 0;
        for (i = 0; i < int(Element_Composition[k].size()); i++) {
            x = x +
                Element_Composition[k][i].prob * Element_Composition[k][i].mw;
            xx = xx + Element_Composition[k][i].prob *
                          Element_Composition[k][i].mw *
                          Element_Composition[k][i].mw;
        }
        v_mw = xx - x * x;
        elements_variance_mass.push_back(v_mw);
    }
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*Computes element monoisotopic weight*/
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double MIDAs::monoisotopic_molecular_mass(
    vector<vector<struct Composition> > &element_composition) {
    int k, i, j;
    double mw = 0.0;
    double min_mw;

    for (k = 0; k < int(element_composition.size()); k++) {
        min_mw = 1e9;
        j = 0;
        for (i = 0; i < int(element_composition[k].size()); i++) {
            // cerr<<Element_Composition[k][i].atoms<<"\t"<<Element_Composition[k][i].mw<<"\t"<<Element_Composition[k][i].prob<<"\t"<<j<<"\t"<<i<<"\t"<<min_mw<<endl;
            if (element_composition[k][i].mw < min_mw) {
                min_mw = element_composition[k][i].mw;
                j = i;
            }
        }
        if (int(element_composition[k].size()) > 0) {
            mw = mw +
                 element_composition[k][j].atoms * element_composition[k][j].mw;
        }
    }
    if (CHARGE_STATE != 0) {
        mw = mw / fabs(CHARGE_STATE);
    }

    return mw;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*Computes theoretical average mass*/
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double MIDAs::average_molecular_mass() {
    int k, i;
    double ave_mw;

    ave_mw = 0.0;

    for (k = 0; k < int(Element_Composition.size()); k++) {
        for (i = 0; i < int(Element_Composition[k].size()); i++) {
            if (CHARGE_STATE != 0) {
                ave_mw = ave_mw + Element_Composition[k][i].atoms *
                                      Element_Composition[k][i].mw *
                                      Element_Composition[k][i].prob /
                                      fabs(CHARGE_STATE);
            } else {
                ave_mw = ave_mw + Element_Composition[k][i].atoms *
                                      Element_Composition[k][i].mw *
                                      Element_Composition[k][i].prob;
            }
        }
    }

    return ave_mw;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*Computes element variances. Returns theoretical variance*/
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double MIDAs::variance_molecular_mass() {
    int k, i;
    double sigma, var_mw, ave_mw;

    sigma = var_mw = ave_mw = 0;
    for (k = 0; k < int(Element_Composition.size()); k++) {
        var_mw = ave_mw = 0;
        if (int(Element_Composition[k].size()) > 0) {
            for (i = 0; i < int(Element_Composition[k].size()); i++) {
                if (CHARGE_STATE != 0) {
                    ave_mw = ave_mw + Element_Composition[k][i].mw *
                                          Element_Composition[k][i].prob /
                                          fabs(CHARGE_STATE);
                } else {
                    ave_mw = ave_mw + Element_Composition[k][i].mw *
                                          Element_Composition[k][i].prob;
                }
            }

            for (i = 0; i < int(Element_Composition[k].size()); i++) {
                if (CHARGE_STATE != 0) {
                    var_mw = var_mw + Element_Composition[k][i].prob *
                                          pow(Element_Composition[k][i].mw /
                                                      fabs(CHARGE_STATE) -
                                                  ave_mw,
                                              2.0);
                } else {
                    var_mw =
                        var_mw +
                        Element_Composition[k][i].prob *
                            pow(Element_Composition[k][i].mw - ave_mw, 2.0);
                }
            }

            var_mw = Element_Composition[k][0].atoms * var_mw;
            sigma = sigma + var_mw;
        }
    }
    return sigma;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*Computes cumulants from probability distribution*/
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void MIDAs::compute_distribution_cumulant_moments(
    vector<struct Isotopic_Distribution> &DISTRIBUTION,
    vector<double> &cumulants, int n) {
    double i;
    double k;
    vector<double> moments;

    // computes distribution moments
    compute_distribution_central_moments(DISTRIBUTION, moments, n);

    // compute distribution cumulant
    for (i = 1; i <= n; i++) {
        if (i == 1) {
            k = moments[0];
            cumulants.push_back(k);
        } else if (i == 2) {
            k = moments[1];
            cumulants.push_back(k);
        } else if (i == 3) {
            k = moments[2];
            cumulants.push_back(k);
        } else if (i == 4) {
            k = moments[3] - 3 * pow(moments[1], 2.0);
            cumulants.push_back(k);
        } else if (i == 5) {
            k = moments[4] - 10 * moments[1] * moments[2];
            cumulants.push_back(k);
        } else if (i == 6) {
            k = moments[5] - 15 * moments[3] * moments[1] -
                10 * pow(moments[2], 2.0) + 30 * pow(moments[1], 3.0);
            cumulants.push_back(k);
        } else if (i == 7) {
            k = moments[6] - 21 * moments[4] * moments[1] -
                35 * moments[3] * moments[2] +
                210 * moments[2] * pow(moments[1], 2.0);
            cumulants.push_back(k);
        } else if (i == 8) {
            k = moments[7] - 28 * moments[5] * moments[1] -
                56 * moments[4] * moments[2] - 35 * pow(moments[3], 2.0) +
                420 * moments[3] * pow(moments[1], 2.0) +
                560 * pow(moments[2], 2.0) * moments[1] -
                630 * pow(moments[1], 4.0);
            cumulants.push_back(k);
        }
    }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*Computes central-moments from probability distribution*/
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void MIDAs::compute_distribution_central_moments(
    vector<struct Isotopic_Distribution> &DISTRIBUTION, vector<double> &moments,
    int n) {
    int i, j;
    double u, ave;

    ave = 0;
    for (i = 0; i < int(DISTRIBUTION.size()); i++) {
        ave = ave + DISTRIBUTION[i].prob * DISTRIBUTION[i].mw;
    }
    moments.push_back(ave);

    for (j = 1; j < n; j++) {
        u = 0;
        for (i = 0; i < int(DISTRIBUTION.size()); i++) {
            u = u + DISTRIBUTION[i].prob * pow(DISTRIBUTION[i].mw - ave, j + 1);
        }
        moments.push_back(u);
    }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*Computes average molecular weight from mass probability distribution*/
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double MIDAs::average_molecular_mass_distribution(
    vector<struct Isotopic_Distribution> &DISTRIBUTION) {
    int i;
    double ave_mw, sum_p;

    // compute average molecular weight
    ave_mw = sum_p = 0.0;
    for (i = 0; i < int(DISTRIBUTION.size()); i++) {
        if (CHARGE_STATE != 0) {
            ave_mw = ave_mw + DISTRIBUTION[i].mw * DISTRIBUTION[i].prob /
                                  fabs(CHARGE_STATE);
        } else {
            ave_mw = ave_mw + DISTRIBUTION[i].mw * DISTRIBUTION[i].prob;
        }
        sum_p = sum_p + DISTRIBUTION[i].prob;
    }
    ave_mw = ave_mw / sum_p;

    return ave_mw;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*Computes variance of the molecular weight from mass probability distribution*/
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double MIDAs::variance_molecular_mass_distribution(
    vector<struct Isotopic_Distribution> &DISTRIBUTION) {
    int i;
    double ave_mw, var_mw, sum_p;

    // computes  molecular weight
    ave_mw = var_mw = sum_p = 0.0;
    for (i = 0; i < int(DISTRIBUTION.size()); i++) {
        if (CHARGE_STATE != 0) {
            ave_mw =
                ave_mw + pow(DISTRIBUTION[i].mw / fabs(CHARGE_STATE), 1.0) *
                             DISTRIBUTION[i].prob;
            var_mw =
                var_mw + pow(DISTRIBUTION[i].mw / fabs(CHARGE_STATE), 2.0) *
                             DISTRIBUTION[i].prob;
        } else {
            ave_mw =
                ave_mw + pow(DISTRIBUTION[i].mw, 1.0) * DISTRIBUTION[i].prob;
            var_mw =
                var_mw + pow(DISTRIBUTION[i].mw, 2.0) * DISTRIBUTION[i].prob;
        }

        sum_p = sum_p + DISTRIBUTION[i].prob;
    }

    ave_mw = ave_mw / sum_p;
    var_mw = var_mw / sum_p;

    var_mw = var_mw - ave_mw * ave_mw;

    // cerr<<"ID:Sum of probabilities = "<<sum_p<<"\t"<<ave_mw<<endl;
    return var_mw;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// returns elemental isotopic composition
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
vector<vector<struct Composition> > MIDAs::ELEMENTAL_ISOTOPIC_COMPOSITION() {
    return Element_Composition;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// returns monoisotopic mass
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double MIDAs::monoisotopic_mass() { return MONOISOTOPIC_MASS; }
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// returns time of FGID
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double MIDAs::fgid_time() { return FGID_TIME; }
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// returns time of CGID
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double MIDAs::cgid_time() { return CGID_TIME; }
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Returns ln(n!)
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double MIDAs::FACTR_LN(int n) {
    static double ntop = 1.0;
    static double a[50003] = {
        0.0, 0.0};  // A static array is automatically initialized to zero.
    int j;

    if (n < 0) {
        cerr << "Negative factorial in routine factln" << endl;
    } else if (n <= 1)
        return 0.0;
    else if (n <= 50000) {
        while (ntop <= n) {
            j = ntop++;
            a[j + 1] = a[j] + log(ntop);
        }
        return a[n];
    }
    return n * log(n) - n + 0.5 * log(6.28318530717959 * n) +
           0.08333333333333 / n - 0.00333333333333 / (n * n * n);
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double MIDAs::GAMMA_LN(double xx) {
    double x, y, tmp, ser;
    static double cof[6] = {76.18009172947146,     -86.50532032941677,
                            24.01409824083091,     -1.231739572450155,
                            0.1208650973866179e-2, -0.5395239384953e-5};
    int j;
    y = x = xx;
    tmp = x + 5.5;
    tmp -= (x + 0.5) * log(tmp);
    ser = 1.000000000190015;
    for (j = 0; j <= 5; j++) ser += cof[j] / ++y;
    return -tmp + log(2.5066282746310005 * ser /
                      x);  // Returns the value ln[ÃƒÂƒ(xx)] for xx > 0.
}
