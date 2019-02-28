#ifndef BPARTICLE_H
#define BPARTICLE_H

#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <string>
#include <vector>
#include <stdint.h>
#include <stdlib.h>
#include "bank.h"

#include "TVector3.h"

class BParticle : public hipo::bank {

	private:

		int pid_order     ;
		int px_order      ;
		int py_order      ;
		int pz_order      ;
		int vx_order      ;
		int vy_order      ;
		int vz_order      ;
		int charge_order  ;
		int beta_order    ;
		int chi2pid_order ;
		int status_order  ;

	public:

		BParticle(){};

		BParticle(const char *bankName, hipo::reader &r) : hipo::bank(bankName,r){
			pid_order     = getEntryOrder("pid"    );
			px_order      = getEntryOrder("px"     );
			py_order      = getEntryOrder("py"     );
			pz_order      = getEntryOrder("pz"     );
			vx_order      = getEntryOrder("vx"     );
			vy_order      = getEntryOrder("vy"     );
			vz_order      = getEntryOrder("vz"     );
			charge_order  = getEntryOrder("charge" );
			beta_order    = getEntryOrder("beta"   );
			chi2pid_order = getEntryOrder("chi2pid");
			status_order  = getEntryOrder("status" );
		}

		~BParticle();

		void  init(const char *bankName, hipo::reader &r);
		int   getPid    (int index) { return getInt  (pid_order    ,index);}
		float getPx     (int index) { return getFloat(px_order     ,index);}
		float getPy     (int index) { return getFloat(py_order     ,index);}
		float getPz     (int index) { return getFloat(pz_order     ,index);}
		float getVx     (int index) { return getFloat(vx_order     ,index);}
                float getVy     (int index) { return getFloat(vy_order     ,index);}
                float getVz     (int index) { return getFloat(vz_order     ,index);}
		int   getCharge (int index) { return getInt  (charge_order ,index);}
		float getBeta   (int index) { return getFloat(beta_order   ,index);}
		float getChi2pid(int index) { return getFloat(chi2pid_order,index);}
		int   getStatus (int index) { return getInt  (status_order ,index);}

		TVector3 getV3v(int index);
		TVector3 getV3P(int index);
};

#endif
