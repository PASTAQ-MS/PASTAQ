/*  $Id: PTMdb.hpp 196490 2010-07-05 21:49:29Z ogurtsov $
 * =============================================================================
 *
 *                            PUBLIC DOMAIN NOTICE
 *               National Center for Biotechnology Information
 *
 *  This software/database is a "United States Government Work" under the
 *  terms of the United States Copyright Act.  It was written as part of
 *  the author's official duties as a United States Government employee and
 *  thus cannot be copyrighted.  This software/database is freely available
 *  to the public for use. The National Library of Medicine and the U.S.
 *  Government have not placed any restriction on its use or reproduction.
 *
 *  Although all reasonable efforts have been taken to ensure the accuracy
 *  and reliability of the software and data, the NLM and the U.S.
 *  Government do not and cannot warrant the performance or results that
 *  may be obtained by using this software or data. The NLM and the U.S.
 *  Government disclaim all warranties, express or implied, including
 *  warranties of performance, merchantability or fitness for any particular
 *  purpose.
 *
 *  Please cite the author in any work or product based on this material.
 *
 * =============================================================================
 *
 * Author:  Aleksey Ogurtsov
 *
 * File Description:
 *   PTMdb.hpp : the database of PTMs.
 *
 */
#include <stdio.h> 
#include <string.h> 
#include <vector>
#include <string>
#include <algorithm>
#include <map>

#if !defined EC_CLASS_TAG
#define EC_CLASS_TAG    1

using namespace std;

enum { xUnknown = 0, xPep = 1, xChem = 2, xDNA = 4, xRNA = 8, xAll = 15 };

class Elem { 
public:
    union { 
        char    name[4]; 
        int     code;
    };
    int count; 
public:
    inline bool operator<( const Elem& r ) const { return code < r.code; };
};
typedef struct {
    char    element[4]; 
    int     atoms;
    double  mw;
    double  prob;
} sComposition;

class ElementComposition {
public:
    typedef struct { const char *id; double mw, p; } Element;
    Element*  elements;
    const char**    pepcmp;
    const char**    modcmp;
    const char**    dnacmp;
    const char**    rnacmp;
    const char*     cm;
    const char*     nt;
    const char*     ct;
    int             xdef;
    class CStr {
    private:
        union {
            char    s[4];
            int     i;
        };
    public:
        inline CStr( const char *_s ) : i(0) { 
          int k;
            if( *_s == ')' && unsigned( _s[1] - '0' ) <= '9' - '0' ) {
                s[0] = ')';
                s[1] = '_';
                return;
            }
            if( unsigned( *_s - 'A' ) > 'Z' - 'A' ) return;
            s[0] = _s[0];
            for( k=1; k<4; k++ ) {
                if( unsigned( _s[k] - '0' ) <= '9' - '0' ) { s[k] = '_'; break; }
                if( unsigned( _s[k] - 'a' ) >  'z' - 'a' )               break;
                s[k] = _s[k];
            }
        };
        inline CStr( const char *_s, int _i ) : i(0) { 
          int k;
            if( *_s == ')' && unsigned( _s[1] - '0' ) <= '9' - '0' ) {
                s[0] = ')';
                s[1] = '_';
                return;
            }
            if( unsigned( *_s - 'A' ) > 'Z' - 'A' ) return;
            s[0] = _s[0];
            for( k=1; k<4; k++ ) {
                if( unsigned( _s[k] - '0' ) <= '9' - '0' ) { s[k] = '_'; break; }
                if( unsigned( _s[k] - 'a' ) >  'z' - 'a' )               break;
                s[k] = _s[k];
            }
            if( _i ) s[k] = '_';
        };
        inline bool operator<( CStr& r ) { return i < r.i; };
        inline operator int( ) { return i; };
        inline operator char *() { return s; };
    };
    map<int, int>   emap;
private:
#if 0
    inline const char* update_elem( const char* e ) {
      int i = 0;
        switch( *e ){
            case 'S': if( e[1] == 'e' ) { i++; e++; }
                      i++;
            case 'O': i++;
            case 'N': i++;
            case 'H': i++;
            case 'C': i++;
        }
        if( i ) {
            count[i-1]++;
            e++;
        } else {
            fprintf( stderr, "Bad element: %s\n", e );
            exit( 1 );
        }
        return e;
    };
    inline const char *update_aa( const char *e ) {

    };
#endif
public:
    inline ElementComposition( const ElementComposition& e ) : elements(e.elements), pepcmp(e.pepcmp), modcmp(e.modcmp), cm(e.cm), nt(e.nt), ct(e.ct) {};
    inline ElementComposition( const char *_cm="C", const char *_nt="H", const char *_ct="OH", int xDef = xPep ) { 
        static Element _elements[] = {
{ "H",	1.007825,	99.9885	},	{ "H",	2.014102,	0.0115	},	{ "H",	3.016049,	0	},	{ "He",	3.016029,	0.000137},
{ "He",	4.002603,	99.999863 },	{ "Li",	6.015122,	7.59	},	{ "Li",	7.016004,	92.41	},	{ "Be",	9.012182,	100	},
{ "B",	10.012937,	19.9	},	{ "B",	11.009305,	80.1	},	{ "C",	12,		98.93	},	{ "C",	13.003355,	1.07	},
{ "C",	14.003242,	0	},	{ "N",	14.003074,	99.632	},	{ "N",	15.000109,	0.368	},	{ "O",	15.994915,	99.757	},
{ "O",	16.999132,	0.038	},	{ "O",	17.99916,	0.205	},	{ "F",	18.998403,	100	},	{ "Ne",	19.99244,	90.48	},
{ "Ne",	20.993847,	0.27	},	{ "Ne",	21.991386,	9.25	},	{ "Na",	22.9877,	100	},	{ "Mg",	23.985042,	78.99	},
{ "Mg",	24.985837,	10	},	{ "Mg",	25.982593,	11.01	},	{ "Al",	26.981538,	100	},	{ "Si",	27.976927,	92.2297	},
{ "Si",	28.976495,	4.6832	},	{ "Si",	29.97377,	3.0872	},	{ "P",	30.973762,	100	},	{ "S",	31.972071,	94.93	},
{ "S",	32.971458,	0.76	},	{ "S",	33.967867,	4.29	},	{ "S",	35.967081,	0.02	},	{ "Cl",	34.968853,	75.78	},
{ "Cl",	36.965903,	24.22	},	{ "Ar",	35.967546,	0.3365	},	{ "Ar",	37.962732,	0.0632	},	{ "Ar",	39.962383,	99.6003	},
{ "K",	38.963707,	93.2581	},	{ "K",	39.963999,	0.0117	},	{ "K",	40.961826,	6.7302	},	{ "Ca",	39.962591,	96.941	},
{ "Ca",	41.958618,	0.647	},	{ "Ca",	42.958767,	0.135	},	{ "Ca",	43.955481,	2.086	},	{ "Ca",	45.953693,	0.004	},
{ "Ca",	47.952534,	0.187	},	{ "Sc",	44.95591,	100	},	{ "Ti",	45.952629,	8.25	},	{ "Ti",	46.951764,	7.44	},
{ "Ti",	47.947947,	73.72	},	{ "Ti",	48.947871,	5.41	},	{ "Ti",	49.944792,	5.18	},	{ "V",	49.947163,	0.25	},
{ "V",	50.943964,	99.75	},	{ "Cr",	49.94605,	4.345	},	{ "Cr",	51.940512,	83.789	},	{ "Cr",	52.940654,	9.501	},
{ "Cr",	53.938885,	2.365	},	{ "Mn",	54.93805,	100	},	{ "Fe",	53.939615,	5.845	},	{ "Fe",	55.934942,	91.754	},
{ "Fe",	56.935399,	2.119	},	{ "Fe",	57.93328,	0.282	},	{ "Co",	58.9332,	100	},	{ "Ni",	57.935348,	68.0769	},
{ "Ni",	59.930791,	26.2231	},	{ "Ni",	60.93106,	1.1399	},	{ "Ni",	61.928349,	3.6345	},	{ "Ni",	63.92797,	0.9256	},
{ "Cu",	62.929601,	69.17	},	{ "Cu",	64.927794,	30.83	},	{ "Zn",	63.929147,	48.63	},	{ "Zn",	65.926037,	27.9	},
{ "Zn",	66.927131,	4.1	},	{ "Zn",	67.924848,	18.75	},	{ "Zn",	69.925325,	0.62	},	{ "Ga",	68.925581,	60.108	},
{ "Ga",	70.924705,	39.892	},	{ "Ge",	69.92425,	20.84	},	{ "Ge",	71.922076,	27.54	},	{ "Ge",	72.923459,	7.73	},
{ "Ge",	73.921178,	36.28	},	{ "Ge",	75.921403,	7.61	},	{ "As",	74.921596,	100	},	{ "Se",	73.922477,	0.89	},
{ "Se",	75.919214,	9.37	},	{ "Se",	76.919915,	7.63	},	{ "Se",	77.91731,	23.77	},	{ "Se",	79.916522,	49.61	},
{ "Se",	81.9167,	8.73	},	{ "Br",	78.918338,	50.69	},	{ "Br",	80.916291,	49.31	},	{ "Kr",	77.920386,	0.35	},
{ "Kr",	79.916378,	2.28	},	{ "Kr",	81.913485,	11.58	},	{ "Kr",	82.914136,	11.49	},	{ "Kr",	83.911507,	57	},
{ "Kr",	85.91061,	17.3	},	{ "Rb",	84.911789,	72.17	},	{ "Rb",	86.909183,	27.83	},	{ "Sr",	83.913425,	0.56	},
{ "Sr",	85.909262,	9.86	},	{ "Sr",	86.908879,	7	},	{ "Sr",	87.905614,	82.58	},	{ "Y",	88.905848,	100	},
{ "Zr",	89.904704,	51.45	},	{ "Zr",	90.905645,	11.22	},	{ "Zr",	91.90504,	17.15	},	{ "Zr",	93.906316,	17.38	},
{ "Zr",	95.908276,	2.8	},	{ "Nb",	92.906378,	100	},	{ "Mo",	91.90681,	14.84	},	{ "Mo",	93.905088,	9.25	},
{ "Mo",	94.905841,	15.92	},	{ "Mo",	95.904679,	16.68	},	{ "Mo",	96.906021,	9.55	},	{ "Mo",	97.905408,	24.13	},
{ "Mo",	99.907477,	9.63	},	{ "Tc",	97.907216,	0	},	{ "Ru",	95.907598,	5.54	},	{ "Ru",	97.905287,	1.87	},
{ "Ru",	98.905939,	12.76	},	{ "Ru",	99.90422,	12.6	},	{ "Ru",	100.905582,	17.06	},	{ "Ru",	101.90435,	31.55	},
{ "Ru",	103.90543,	18.62	},	{ "Rh",	102.905504,	100	},	{ "Pd",	101.905608,	1.02	},	{ "Pd",	103.904035,	11.14	},
{ "Pd",	104.905084,	22.33	},	{ "Pd",	105.903483,	27.33	},	{ "Pd",	107.903894,	26.46	},	{ "Pd",	109.905152,	11.72	},
{ "Ag",	106.905093,	51.839	},	{ "Ag",	108.904756,	48.161	},	{ "Cd",	105.906458,	1.25	},	{ "Cd",	107.904183,	0.89	},
{ "Cd",	109.903006,	12.49	},	{ "Cd",	110.904182,	12.8	},	{ "Cd",	111.902757,	24.13	},	{ "Cd",	112.904401,	12.22	},
{ "Cd",	113.903358,	28.73	},	{ "Cd",	115.904755,	7.49	},	{ "In",	112.904061,	4.29	},	{ "In",	114.903878,	95.71	},
{ "Sn",	111.904821,	0.97	},	{ "Sn",	113.902782,	0.66	},	{ "Sn",	114.903346,	0.34	},	{ "Sn",	115.901744,	14.54	},
{ "Sn",	116.902954,	7.68	},	{ "Sn",	117.901606,	24.22	},	{ "Sn",	118.903309,	8.59	},	{ "Sn",	119.902197,	32.58	},
{ "Sn",	121.90344,	4.63	},	{ "Sn",	123.905275,	5.79	},	{ "Sb",	120.903818,	57.21	},	{ "Sb",	122.904216,	42.79	},
{ "Te",	119.90402,	0.09	},	{ "Te",	121.903047,	2.55	},	{ "Te",	122.904273,	0.89	},	{ "Te",	123.902819,	4.74	},
{ "Te",	124.904425,	7.07	},	{ "Te",	125.903306,	18.84	},	{ "Te",	127.904461,	31.74	},	{ "Te",	129.906223,	34.08	},
{ "I",	126.904468,	100	},	{ "Xe",	123.905896,	0.09	},	{ "Xe",	125.904269,	0.09	},	{ "Xe",	127.90353,	1.92	},
{ "Xe",	128.904779,	26.44	},	{ "Xe",	129.903508,	4.08	},	{ "Xe",	130.905082,	21.18	},	{ "Xe",	131.904154,	26.89	},
{ "Xe",	133.905395,	10.44	},	{ "Xe",	135.90722,	8.87	},	{ "Cs",	132.905447,	100	},	{ "Ba",	129.90631,	0.106	},
{ "Ba",	131.905056,	0.101	},	{ "Ba",	133.904503,	2.417	},	{ "Ba",	134.905683,	6.592	},	{ "Ba",	135.90457,	7.854	},
{ "Ba",	136.905821,	11.232	},	{ "Ba",	137.905241,	71.698	},	{ "La",	137.907107,	0.09	},	{ "La",	138.906348,	99.91	},
{ "Ce",	135.907144,	0.185	},	{ "Ce",	137.905986,	0.251	},	{ "Ce",	139.905434,	88.45	},	{ "Ce",	141.90924,	11.114	},
{ "Pr",	140.907648,	100	},	{ "Nd",	141.907719,	27.2	},	{ "Nd",	142.90981,	12.2	},	{ "Nd",	143.910083,	23.8	},
{ "Nd",	144.912569,	8.3	},	{ "Nd",	145.913112,	17.2	},	{ "Nd",	147.916889,	5.7	},	{ "Nd",	149.920887,	5.6	},
{ "Pm",	144.912744,	0	},	{ "Sm",	143.911995,	3.07	},	{ "Sm",	146.914893,	14.99	},	{ "Sm",	147.914818,	11.24	},
{ "Sm",	148.91718,	13.82	},	{ "Sm",	149.917271,	7.38	},	{ "Sm",	151.919728,	26.75	},	{ "Sm",	153.922205,	22.75	},
{ "Eu",	150.919846,	47.81	},	{ "Eu",	152.921226,	52.19	},	{ "Gd",	151.919788,	0.2	},	{ "Gd",	153.920862,	2.18	},
{ "Gd",	154.922619,	14.8	},	{ "Gd",	155.92212,	20.47	},	{ "Gd",	156.923957,	15.65	},	{ "Gd",	157.924101,	24.84	},
{ "Gd",	159.927051,	21.86	},	{ "Tb",	158.925343,	100	},	{ "Dy",	155.924278,	0.06	},	{ "Dy",	157.924405,	0.1	},
{ "Dy",	159.925194,	2.34	},	{ "Dy",	160.92693,	18.91	},	{ "Dy",	161.926795,	25.51	},	{ "Dy",	162.928728,	24.9	},
{ "Dy",	163.929171,	28.18	},	{ "Ho",	164.930319,	100	},	{ "Er",	161.928775,	0.14	},	{ "Er",	163.929197,	1.61	},
{ "Er",	165.93029,	33.61	},	{ "Er",	166.932045,	22.93	},	{ "Er",	167.932368,	26.78	},	{ "Er",	169.93546,	14.93	},
{ "Tm",	168.934211,	100	},	{ "Yb",	167.933894,	0.13	},	{ "Yb",	169.934759,	3.04	},	{ "Yb",	170.936322,	14.28	},
{ "Yb",	171.936378,	21.83	},	{ "Yb",	172.938207,	16.13	},	{ "Yb",	173.938858,	31.83	},	{ "Yb",	175.942568,	12.76	},
{ "Lu",	174.940768,	97.41	},	{ "Lu",	175.942682,	2.59	},	{ "Hf",	173.94004,	0.16	},	{ "Hf",	175.941402,	5.26	},
{ "Hf",	176.94322,	18.6	},	{ "Hf",	177.943698,	27.28	},	{ "Hf",	178.945815,	13.62	},	{ "Hf",	179.946549,	35.08	},
{ "Ta",	179.947466,	0.012	},	{ "Ta",	180.947996,	99.988	},	{ "W",	179.946706,	0.12	},	{ "W",	181.948206,	26.5	},
{ "W",	182.950224,	14.31	},	{ "W",	183.950933,	30.64	},	{ "W",	185.954362,	28.43	},	{ "Re",	184.952956,	37.4	},
{ "Re",	186.955751,	62.6	},	{ "Os",	183.952491,	0.02	},	{ "Os",	185.953838,	1.59	},	{ "Os",	186.955748,	1.96	},
{ "Os",	187.955836,	13.24	},	{ "Os",	188.958145,	16.15	},	{ "Os",	189.958445,	26.26	},	{ "Os",	191.961479,	40.78	},
{ "Ir",	190.960591,	37.3	},	{ "Ir",	192.962924,	62.7	},	{ "Pt",	189.95993,	0.014	},	{ "Pt",	191.961035,	0.782	},
{ "Pt",	193.962664,	32.967	},	{ "Pt",	194.964774,	33.832	},	{ "Pt",	195.964935,	25.242	},	{ "Pt",	197.967876,	7.163	},
{ "Au",	196.966552,	100	},	{ "Hg",	195.965815,	0.15	},	{ "Hg",	197.966752,	9.97	},	{ "Hg",	198.968262,	16.87	},
{ "Hg",	199.968309,	23.1	},	{ "Hg",	200.970285,	13.18	},	{ "Hg",	201.970626,	29.86	},	{ "Hg",	203.973476,	6.87	},
{ "Tl",	202.972329,	29.524	},	{ "Tl",	204.974412,	70.476	},	{ "Pb",	203.973029,	1.4	},	{ "Pb",	205.974449,	24.1	},
{ "Pb",	206.975881,	22.1	},	{ "Pb",	207.976636,	52.4	},	{ "Bi",	208.980383,	100	},	{ "Po",	208.982416,	0	},
{ "At",	209.987131,	0	},	{ "Rn",	222.01757,	0	},	{ "Fr",	223.019731,	0	},	{ "Ra",	226.025403,	0	},
{ "Ac",	227.027747,	0	},	{ "Th",	232.03805,	100	},	{ "Pa",	231.035879,	100	},	{ "U",	234.040946,	0.0055	},
{ "U",	235.043923,	0.72	},	{ "U",	238.050783,	99.2745	},	{ "Np",	237.048167,	0	},	{ "Pu",	244.064198,	0	},
{ "Am",	243.061373,	0	},	{ "Cm",	247.070347,	0	},	{ "Bk",	247.070299,	0	},	{ "Cf",	251.07958,	0	},
{ "Es",	252.082972,	0	},	{ "Fm",	257.095099,	0	},	{ "Md",	258.098425,	0	},	{ "No",	259.101024,	0	},
{ "Lr",	262.109692,	0	},	{ "Rf",	263.118313,	0	},	{ "Db",	262.011437,	0	},	{ "Sg",	266.012238,	0	},
{ "Bh",	264.012496,	0	},	{ "Hs",	269.001341,	0	},	{ "Mt",	268.001388,	0	},	{ "Uun",272.001463,	0	},
{ "Uuu",272.001535,	0	},	{ "Uub",277,		0	},	{ "Uuq",	289,	0	},	{ "Uuh",	289,	0	},
{ "Uuo",	293,	0 	},  { "",   0,          0   }, };
        static const char*  _pepcmp[] = {
    "Trp", "C11H10N2O", "Pyl", "C12H19N3O2", "Gly", "C2H3NO",  "Ala", "C3H5NO",   "Cys", "C3H5NOS",  "Sec", "C3H5NOSe", "Ser", "C3H5NO2",
    "Asp", "C4H5NO3",   "Asn", "C4H6N2O2",   "Thr", "C4H7NO2", "Gln", "C5H8N2O2", "Tyr", "C9H9NO2",  "Val", "C5H9NO",   "Met", "C5H9NOS",
    "Pro", "C5H7NO",    "Glu", "C5H7NO3",    "Ile", "C6H11NO", "Leu", "C6H11NO",  "Lys", "C6H12N2O", "Arg", "C6H12N4O", "His", "C6H7N3O",
    "Phe", "C9H9NO",    "Ac", "CH3CO", 
    "W", "C11H10N2O", "O", "C12H19N3O2", "G", "C2H3NO",  "A", "C3H5NO",   "C", "C3H5NOS",  "U", "C3H5NOSe", "S", "C3H5NO2",
    "D", "C4H5NO3",   "N", "C4H6N2O2",   "T", "C4H7NO2", "Q", "C5H8N2O2", "Y", "C9H9NO2",  "V", "C5H9NO",   "M", "C5H9NOS",
    "P", "C5H7NO",    "E", "C5H7NO3",    "I", "C6H11NO", "L", "C6H11NO",  "K", "C6H12N2O", "R", "C6H12N4O", "H", "C6H7N3O",
    "F", "C9H9NO",   
    "", "" };
        static const char*  _dnacmp[] = { "A", "C10N5O5H12P", "C", "C8N3O6H12P", "G", "C10N5O6N12P", "T", "C10N2O7H13P", "", "" };
        static const char*  _rnacmp[] = { "A", "C10N5O6H12P", "C", "C8N3O7H12P", "G", "C10N5O7N12P", "U", "C7N2O8H11P", "T", "C7N2O8H11P", "", "" };
        static const char*  _modcmp[] = { "O", "O", "M", "CH2", "C", "OCH2", "S", "SO3", "P", "OPO3H", "H", "O", "", "" }; 
               const char*  _cysmod[] = { "C00", "NC3H5OS", "C31", "NC5H7O3S", "C32", "N2C5H8O2S", "C33", "N2C10H11OS", "", "" };  
        elements = _elements;
        pepcmp   = _pepcmp;
        modcmp   = _modcmp;
        dnacmp   = _dnacmp;
        rnacmp   = _rnacmp;
        xdef     =  xDef;
        ct       = _ct;
        nt       = _nt;
        cm       = _cm;
        const char **s;
        if( cm[0] == 'C' && digit( cm[1] ) ) for( s=_cysmod; !!**s; s+=2 ) if( ex( *s, cm ) ) {
            for( const char** c=_pepcmp; !!**c; c+=2 ) if( ex( *c, "C" ) || ex( *c, "Cys" ) ) c[1] = s[1];
            break;
        }
        for( Element* e=_elements; !!e->id[0]; e++ ) {
            emap[CStr(e->id)] = xChem;
            emap[CStr(e->id,1)] = xChem;
        }
        emap[CStr(")0")] |= xChem;
        for( s=_pepcmp; **s; s+=2 ) emap[CStr(*s)] |= xPep;
        for( s=_dnacmp; **s; s+=2 ) emap[CStr(*s)] |= xDNA;
        for( s=_rnacmp; **s; s+=2 ) emap[CStr(*s)] |= xRNA;
    };
    inline bool uppercase( char c ) { return (unsigned)( c - 'A' ) <= 'Z' - 'A'; };
    inline bool lowercase( char c ) { return (unsigned)( c - 'a' ) <= 'z' - 'a'; };
    inline bool digit( char c ) { return (unsigned)( c - '0' ) <= '9' - '0'; };
    inline bool in( const char *s, char c ) {
        while( !!*s ) if( *s++ == c ) return true;
        return false;
    }
    inline int typeone( const char *p ) {
      int t = xUnknown;
        if( in( "ADEGLMQRT", *p ) && !lowercase( p[1] ) ) t |= xPep;
        if( digit( p[1] ) ) t |= xChem;
        if( in( "CHMOPS", *p ) && p[1] == '(' ) t |= xPep;
        if( *p == 'U' && p[1] == 'u' && p[2] ) t |= xChem;
        else if( lowercase( p[1] ) && lowercase( p[2] ) ) t |= xPep;
        return t;
    };
    inline int type_imp( const char *p ) {
      int t = xAll;
        while( *p && *p != '-' && *p != '=' ) {
            if( uppercase( *p ) ) t &= emap[CStr(p)];
            p++;
        }
        return t;
    };
    inline int type( const char *p ) { return type_imp( p ); };
    inline int nttype( const char *p ) {
      int t = xAll;
        while( *p && *p != '-' && *p != '=' ) {
            if( p[0] == 'A' && p[1] == 'c' && !lowercase( p[2] ) ) t &= xPep; else t &= emap[CStr(p)];
            p++;
        }
        return t & xChem ? xChem : t;
    };
    inline int cttype( const char *p ) { 
        int t = type_imp( p );
        return t & xChem ? xChem : t;
    };
    inline bool ex( const char *p, const char *s ) {
        while( !!*s && !!*p ) if( (*s++^*p++) != 0 ) return false;
        return *s == *p;
    };
    inline bool eq( const char *p, const char *s ) {
        while( !!*s && !!*p ) if( (*s++^*p++) != 0 ) return false;
        return *s == 0;
    };
    inline bool update( vector<Elem>& v, const char *p, int t ) {
      size_t i, j, c; Elem a;
//        if( (t&xdef) != 0 ) return false;
        if( t == xUnknown ) t = xPep;
        string s = "";
        if( t == xPep ) {
            for( ; *p&&*p!='-'&&*p!='='; p++ ) {
                if( uppercase( *p ) ) {
                    if( p[1] == '(' ) {
                        for( i=0; modcmp[i][0]!=0; i+=2 ) if( eq( p, modcmp[i] ) ) { s += modcmp[i + 1]; break; }
                        if( modcmp[i][0]==0 ) return false;
                    } else {
                        for( i=0; pepcmp[i][0]!=0; i+=2 ) if( eq( p, pepcmp[i] ) ) { s += pepcmp[i + 1]; break; }
                        if( pepcmp[i][0]==0 ) return false;
                    }
                }
            }
        } else if( t == xDNA ) {
            for( ; *p&&*p!='-'&&*p!='='; p++ ) for( i=0; dnacmp[i][0]!=0; i+=2 ) if( eq( p, dnacmp[i] ) ) { s += dnacmp[i + 1]; break; }
        } else if( t == xRNA ) {
            for( ; *p&&*p!='-'&&*p!='='; p++ ) for( i=0; rnacmp[i][0]!=0; i+=2 ) if( eq( p, rnacmp[i] ) ) { s += rnacmp[i + 1]; break; }
        } else if( t == xChem ) {
            for( ; *p&&*p!='-'&&*p!='='; p++ ) s.push_back( *p );
        } else {
          vector<const char *> v;
            if( t&xChem ) v.push_back( "chemical formula" );
            if( t&xPep  ) v.push_back( "amino acid sequence" );
            if( t&xDNA  ) v.push_back( "DNA" );
            if( t&xRNA  ) v.push_back( "RNA" );
            if( v.size() == 0 ) {
                fprintf( stderr, "Unknown type.\n" );
            } else {
                fprintf( stderr, "Unknown type. Recognized as " );
                for( i=0; i<v.size(); i++ ) {
                    if( i == 0            ) fprintf( stderr,     "%s",  v[i] ); else
                    if( i == v.size() - 1 ) fprintf( stderr, " or %s.", v[i] ); else
                                            fprintf( stderr,   ", %s",  v[i] ); 
                }
                v.clear();
                fprintf( stderr, " Requested type is " );
                if( xdef == xAll ) v.push_back( "any" ); else {
                    if( xdef&xChem ) v.push_back( "chemical formula" );
                    if( xdef&xPep  ) v.push_back( "amino acid sequence" );
                    if( xdef&xDNA  ) v.push_back( "DNA" );
                    if( xdef&xRNA  ) v.push_back( "RNA" );
                }
                for( i=0; i<v.size(); i++ ) {
                    if( i == 0            ) fprintf( stderr,     "%s",  v[i] ); else
                    if( i == v.size() - 1 ) fprintf( stderr, " or %s.", v[i] ); else
                                            fprintf( stderr,   ", %s",  v[i] ); 
                }
                fprintf( stderr, ".\n" );
            }
            exit( 2 );
        }
        vector<Elem> ve;
        vector<size_t> vi;
        for( i=0; s[i]&&s[i]!='-'&&s[i]!='='; ) {
            if( s[i] == '(' ) {
                vi.push_back( ve.size() );
                i++;
            } else if( s[i] == ')' ) {
                if( vi.size() == 0 ) return false;
                if( (c = atoi( &s[i + 1] )) == 0 ) c = 1;
                if( c > 1 ) for( j=vi[vi.size()-1]; j<ve.size(); j++ ) ve[j].count *= c;
                vi.pop_back();
                for( i++; s[i]; i++ ) if( !digit( s[i] ) ) break;
            } else if( uppercase( s[i] ) ) {
                for( a.code=s[i++],c=8; lowercase(s[i]); i++,c+=8 ) a.code += s[i]<<c;
                if( (a.count = atoi( &s[i] )) == 0 ) a.count = 1;
                ve.push_back( a );
            } else {
                i++;
            }
        }
        for( i=0; i<v.size(); i++ ) {
            a.code  = v[i].code;
            a.count = v[i].count;
            ve.push_back( a );
        }
        sort( ve.begin(), ve.end() );
        v.clear();
        for( i=j=0; j<ve.size(); j++ ) {
            if( i < v.size() && ve[j].code == v[i].code ) {
                v[i].count += ve[j].count;
            } else {
                a.code  = ve[j].code;
                a.count = ve[j].count;
                i = v.size();
                v.push_back( a );
            }
        }
        return true;
    };
    template<class User> inline bool CalcComposition( const char *pep, vector< vector<User> >& v ){
      vector<const char *> stp;
      vector<int         > stt;
      size_t               i, j;
      int                  tc, tn;
      vector<Elem> comp;
        v.clear();
        tc = tn = 0;
        for( const char *s=pep; *s; ) {
            stp.push_back( s );
            while( *s && *s != '-' && *s != '=' ) s++;
            if( *s ) s++;
        }
        i = type( pep );
        if( stp.size() == 1 && !!(i&(xPep|xDNA|xRNA)) ) {
            stp.clear();
            if( !!nt && !!*nt ) { stp.push_back( nt  ); tn = 1; }
            stp.push_back( pep );
            if( !!ct && !!*ct ) { stp.push_back( ct  ); tc = 1; }
        }
        for( j=xUnknown,i=0; i<stp.size(); i++ ) {
          int t;
            j |= t = type( stp[i] );
            stt.push_back( t );
        }
        if( (j&(xPep|xDNA|xRNA)) != 0 ) {
            if( tn != 0 ) { i = 0;              if( (stt[i]&xChem) == xChem ) stt[i] = nttype( stp[i] ); }
            for( i=tn; i<stt.size()-tc; i++ ) if( (stt[i]&xdef) != 0 ) stt[i] = j&xdef;
            if( tc != 0 ) { i = stt.size() - 1; if( (stt[i]&xChem) == xChem ) stt[i] = cttype( stp[i] ); }
        } else {
            for( i=0; i<stt.size()  ; i++ ) if( (stt[i]&xdef) != 0 ) stt[i] = j&xdef;
//          for( i=0; i<stt.size()  ; i++ ) stt[i] = xChem;
        }
///        if( stp.size() == 3 ) {
///            if( !update( comp, stp[0], nttype( stp[0] ) ) ) return false;
///            if( !update( comp, stp[1],   type( stp[1] ) ) ) return false;
///            if( !update( comp, stp[2], cttype( stp[2] ) ) ) return false;
///        } else if( stp.size() > 3 ) {
///            for( i=0; i<stp.size(); i++ ) {
///                if( i == 0              && !update( comp, stp[i], nttype( stp[i] ) ) ) return false; else 
///                if( i == stp.size() - 1 && !update( comp, stp[i], cttype( stp[i] ) ) ) return false; else
///                if(                        !update( comp, stp[i],   type( stp[i] ) ) ) return false; 
///            }
///        } else {
            for( i=0; i<stp.size(); i++ ) if( !update( comp, stp[i], stt[i] ) ) return false;
///        }
        for( i=0; i<comp.size(); i++ ) {
            int c = comp[i].code; 
            char* s;
            for( s=comp[i].name; c!=0; c>>=8 ) *s++ = c&255;
            *s = 0;
            v.push_back( vector<User> () );
            vector<User>& v2 = v[v.size()-1];
            for( const Element *e=elements; e->id[0]!=0 ; e++ ) {
                if( ex( comp[i].name, e->id ) && e->p != 0 ) {
                    int c = comp[i].code; 
                    User u;
                    for( s=u.element; c!=0; c>>=8 ) *s++ = c&255;
                    *s = 0;
                    u.atoms = comp[i].count;
                    u.mw    = e->mw;
                    u.prob  = e->p;
                    v2.push_back( u );
                }
            }
        }
        return true;
    };
    inline bool UpdateElement( const char *name, double mw, double prob ) {
        for( Element *e=elements; e->id[0]!=0 ; e++ ) if( ex( name, e->id ) && e->mw >= mw - 0.4 && e->mw <= mw + 0.4 ) {
            e->mw = mw;
            e->p  = prob;
            return true;
        }
        return false;
    };
    inline bool UpdateElementMw( const char *name, double mw ) {
        for( Element *e=elements; e->id[0]!=0 ; e++ ) if( ex( name, e->id ) && e->mw >= mw - 0.4 && e->mw <= mw + 0.4 ) {
            e->mw = mw;
            return true;
        }
        return false;
    };
    inline bool UpdateElementProb( const char *name, double mw, double prob ) {
        for( Element *e=elements; e->id[0]!=0 ; e++ ) if( ex( name, e->id ) && e->mw >= mw - 0.4 && e->mw <= mw + 0.4 ) {
            e->p = prob;
            return true;
        }
        return false;
    };
    template<class User> inline bool UpdateElement    ( User& u ) { return UpdateElement    ( u.element, u.mw, u.prob ); };
    template<class User> inline bool UpdateElementMw  ( User& u ) { return UpdateElementMw  ( u.element, u.mw         ); };
    template<class User> inline bool UpdateElementProb( User& u ) { return UpdateElementProb( u.element, u.mw, u.prob ); };
    inline ~ElementComposition() {};
};

#endif
