#ifndef SPLINTER_MATLAB_H
#define SPLINTER_MATLAB_H

#define API __declspec(dllexport)

#ifdef __cplusplus
	extern "C"
	{
#endif
		void API init();
		double API eval(double x);
#ifdef __cplusplus
	}
#endif


#endif // SPLINTER_MATLAB_H