#ifndef _qk_lib_functions_H
#define _qk_lib_functions_H

#define SAFE_DELETE(p) if(p){delete p; p=NULL;}
#define SAFE_DELETE_ARRAY(p) if(p){delete []p; p=NULL;}

#endif // _qk_lib_functions_H
