#ifndef BBAND_H
#define BBAND_H

#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <string>
#include <vector>
#include <stdint.h>
#include <stdlib.h>
#include "bank.h"

#include "TVector3.h"

class BBand : public hipo::bank {

	private:

		int id_order           ;
		int sector_order       ;
		int layer_order        ;
		int component_order    ;
		int meantimeTdc_order  ;
		int meantimeFadc_order ;
		int difftimeTdc_order  ;
		int difftimeFadc_order ;
		int adcLcorr_order     ;
		int adcRcorr_order     ;
		int tFadcLcorr_order   ;
		int tFadcRcorr_order   ;
		int tTdcLcorr_order    ;
		int tTdcRcorr_order    ;
		int x_order            ;
		int y_order            ;
		int z_order            ;
		int ux_order           ;
		int uy_order           ;
		int uz_order           ;

	public:

		BBand(){};

		BBand(hipo::schema __schema);

		~BBand();


		int   getId           (int index) { return getInt   (id_order           ,index);}
		int   getSector       (int index) { return getInt   (sector_order       ,index);}
		int   getLayer        (int index) { return getInt   (layer_order        ,index);}
		int   getComponent    (int index) { return getInt   (component_order    ,index);}
		float getMeantimeTdc  (int index) { return getFloat (meantimeTdc_order  ,index);}
		float getMeantimeFadc (int index) { return getFloat (meantimeFadc_order ,index);}
		float getDifftimeTdc  (int index) { return getFloat (difftimeTdc_order  ,index);}
		float getDifftimeFadc (int index) { return getFloat (difftimeFadc_order ,index);}
		float getAdcLcorr     (int index) { return getFloat (adcLcorr_order     ,index);}
		float getAdcRcorr     (int index) { return getFloat (adcRcorr_order     ,index);}
		float getTfadcLcorr   (int index) { return getFloat (tFadcLcorr_order   ,index);}
		float getTfadcRcorr   (int index) { return getFloat (tFadcRcorr_order   ,index);}
		float getTtdcLcorr    (int index) { return getFloat (tTdcLcorr_order    ,index);}
		float getTtdcRcorr    (int index) { return getFloat (tTdcRcorr_order    ,index);}
		float getX            (int index) { return getFloat (x_order            ,index);}
		float getY            (int index) { return getFloat (y_order            ,index);}
		float getZ            (int index) { return getFloat (z_order            ,index);}
		float getUx           (int index) { return getFloat (ux_order           ,index);}
		float getUy           (int index) { return getFloat (uy_order           ,index);}
		float getUz           (int index) { return getFloat (uz_order           ,index);}

		int   getBarKey(int index);
};

#endif
