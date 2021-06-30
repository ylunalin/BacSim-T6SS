#ifndef MUTEX_HH

#ifdef _OPENMP
# include <omp.h>
struct mutex {
	mutex() {omp_init_lock(&l);}
	~mutex() {omp_destroy_lock(&l);}
	inline void lock() {omp_set_lock(&l);}
	inline void unlock() {omp_unset_lock(&l);}
	mutex(const mutex& ) {omp_init_lock(&l);}
	mutex& operator= (const mutex& ) {return *this;}
	omp_lock_t l;
};
#else
struct mutex {
	inline void lock() {}
	inline void unlock() {}
};
#endif

#endif
