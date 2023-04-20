//------------------------------------------
//@author KatsuyaKikuchi
//@brief Implementation of the InputKey.h
//------------------------------------------

#include"InputKey.h"

namespace KeyBoard{
	int Key[256] = { 0 };

	int gpUpdateKey(){
		char tmpKey[256];
		GetHitKeyStateAll(tmpKey);
		for (int i = 0; i < 256; i++){
			if (tmpKey[i] != 0){
				Key[i]++;
			}
			else {
				Key[i] = 0;
			}
		}
		return 0;
	}
}
