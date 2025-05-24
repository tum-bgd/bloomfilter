#include <pybind11/pybind11.h>
#include<pybind11/numpy.h>
#include<pybind11/stl.h>

#include <boost/config.hpp>
#include <boost/mpl/assert.hpp>
#include <boost/mpl/int.hpp>

#include <boost/geometry.hpp>
#include <boost/geometry/core/access.hpp>
#include <boost/geometry/core/assert.hpp>
#include <boost/geometry/core/coordinate_type.hpp>
#include <boost/geometry/core/coordinate_system.hpp>
#include <boost/geometry/core/coordinate_dimension.hpp>

#include<iostream>
#include<array>
#include "globimap.hpp"



int add(int i, int j) {
    return i + j;
}
namespace py = pybind11;
namespace bg = boost::geometry;

template<typename func>
void map_matrix(py::array_t<double> mat, func f){
    auto ref = mat.unchecked<2>();
    auto shapeX = ref.shape(0);
    auto shapeY = ref.shape(1);
//    auto dataPtr = (double*) mat.request().ptr;

    for (int x=0; x < shapeX; x++)
        for(int y=0; y < shapeY; y++)
        {
            
            f(x,y,ref(x,y));
        }
}


typedef GloBiMap<bool> globimap_t;
class globimap3d_t :public GloBiMap<bool> {int bogus;};


/*New Bloom Filter Interface


*/
class BloomFilterFacade{
   private:

    inline void hash(const void *data, const size_t len, uint64_t *v1, uint64_t *v2)
    {
	uint64_t hash[2];
	trajcomp::murmur::MurmurHash3_x86_128 ( data, len, *v1, (void*) hash ); // len*sizeof(uint64_t) 
	*v1 = hash[0];
	*v2 = hash[1];
    }


    


   
      std::vector<unsigned char> row_as_array(const py::array &mat, int row ){
          std::vector<unsigned char> hash_data;
	  auto data = mat.data();
	  hash_data.clear();
	  for (int col=0; col < mat.shape(1); col++){
		auto offset = mat.offset_at(row,col);
		const char *p = static_cast<const char *>(data+offset);
		for (size_t i=0; i < mat.itemsize(); i++)
		   hash_data.push_back(*p++);
	    }
	 return hash_data;
      }

   public:
     size_t k,m;
     std::vector<bool> filter;

     std::vector<unsigned char> debug_getrow(const py::array &mat, int row_index){
        auto row = row_as_array(mat,row_index);
	return row;
    }

    std::string debug_gethex(const py::array &mat, int row_index){
        auto row = row_as_array(mat,row_index);
	std::stringstream ss;
	ss << std::hex;
	for (auto c:row)
	  ss << (int) c; 
	return ss.str();
    }
   
    void configure (size_t _k, size_t _m){
          std::cout << "configuring " << _k <<"," << _m << std::endl;
	  k = _k; m=_m;
	  filter.resize(m);
	  for (size_t i=0; i < filter.size(); i++)
	      filter[i] = 0;
      }
    void insert(py::array mat, int axis=0){
	// axis not yet supported, maybe never
	std::vector<unsigned char> hash_data;
        auto data = mat.data();
	for (int item =0; item < mat.shape(0); item ++){
	   hash_data = row_as_array(mat,item);
	    uint64_t h1=8589845122,h2=8465418721;	    
	    hash (&hash_data[0], hash_data.size() ,&h1, &h2);
	    for (size_t i=0; i < static_cast<size_t>(k); i++)
	    {
		   uint64_t v = (h1 + (i+1)*h2)  % m;
		    filter[v] = 1;
	    }
	   }
	
    }


    std::vector<bool> test(py::array mat, int axis=0){
        std::vector<bool> ret(mat.shape(0));
	
	std::vector<unsigned char> hash_data;
        auto data = mat.data();
	for (int item =0; item < mat.shape(0); item ++){
	   ret [item] = true;
	   
	   hash_data = row_as_array(mat,item);
	    uint64_t h1=8589845122,h2=8465418721;	    
	    hash (&hash_data[0], hash_data.size() ,&h1, &h2);
	    for (size_t i=0; i < static_cast<size_t>(k); i++)
	    {
		   uint64_t v = (h1 + (i+1)*h2)  % m;
		   //std::cout << "test hash id: " <<i<<"="<< v << std::endl;
		    if (!filter[v]){
		      ret[item] = false;
		      break;
		   }
	    }
	   }
	 return ret;

    }


    void tobuffer(std::vector<char> &buf)
    {
	buf.clear();
	char ch=0;
	for (size_t i=0; i<filter.size(); i++)
	{
	    auto bit = i %8;
	    ch |= (filter[i] << bit);
	    if (bit == 7){
		// full byte has been built.
		buf.push_back(ch);
		ch = 0;
	    }
	}
	if (filter.size() % 8 != 0)  // we have collected beyond the end.
	  buf.push_back(ch);
    }
    
    void frombuffer(const std::vector<char> &buf, size_t _k, size_t n)
    {
	this->k = _k;
	this->m = n; 
	filter.resize(n);
	size_t k=0;
	for (size_t i=0; i < buf.size(); i++)
	{
	    char ch = buf[i];
	    for(size_t j=0; j < 8; j++)
	      if (k < n)
	        filter[k++] = ch & (1<<j);
	
	}
    }


};




PYBIND11_MODULE(bloomfilter, m) {
    /**
     * NOT FOR PRODUCTIVITY - DEVELOPMENT ONLY
     * */

    m.doc() = R"pbdoc(
        Pybind11 example plugin
        -----------------------

        .. currentmodule:: python_example

        .. autosummary::
           :toctree: _generate

           add
           subtract
    )pbdoc";
    /*Bloom Filter*/
    py::class_<BloomFilterFacade>(m, "bloomfilter",py::module_local())
            .def(py::init<>())
            .def("configure", +[](BloomFilterFacade &self, size_t k, size_t m)  {  self.configure(k,m);})
            .def("configure_fp", +[](BloomFilterFacade &self, size_t n, doube fpr)  {
	        //
		double m = -n*log(fpr) / (log(2)*log(2));
	    self.configure(k,m);
	    })
            .def("insert", +[](BloomFilterFacade &self, py::array &a)  {  self.insert(a);})
	    .def("test", +[](BloomFilterFacade &self, py::array &a)  {  return self.test(a);})
	    .def("insert_along", +[](BloomFilterFacade &self, py::array &a, int axis)  {  self.insert(a,axis);})
	    .def("dbg_getrow",   +[](BloomFilterFacade &self, py::array &a, int row){ return self.debug_getrow(a,row);})
	    .def("dbg_gethex",   +[](BloomFilterFacade &self, py::array &a, int row){ return self.debug_gethex(a,row);})
	    .def("to_bytes",    +[](BloomFilterFacade &self) {std::vector<char> out; self.tobuffer(out);return py::bytes(out.data(),out.size());})
	    .def("from_bytes", +[](BloomFilterFacade &self, py::bytes bytedata, int k, int m){
		 py::buffer_info info(py::buffer(bytedata).request());
                 const char *data = reinterpret_cast<const char *>(info.ptr);
                 size_t length = static_cast<size_t>(info.size);
                 self.frombuffer(std::vector<char>(data, data + length),k,m);
		})
	    
	    ;
//    m.def("test", &test,R"pbdoc(
//        A test for numpy
//
//        Some other explanation about the add function.
//    )pbdoc");
    m.def("add", &add, R"pbdoc(
        Add two numbers

        Some other explanation about the add function.
    )pbdoc");

    m.def("subtract", [](int i, int j) { return i - j; }, R"pbdoc(
        Subtract two numbers

        Some other explanation about the subtract function.
    )pbdoc");

    /**
    * END  - NOT FOR PRODUCTIVITY - DEVELOPMENT ONLY
    * */

#ifdef VERSION_INFO
        m.attr("__version__") = VERSION_INFO;
#else
        m.attr("__version__") = "dev";
#endif

    py::class_<globimap_t>(m, "globimap",py::module_local())
            .def(py::init<>())
            .def("configure", +[](globimap_t &self, size_t k, size_t m)  {  self.configure(k,m);})
            .def("summary", +[](globimap_t &self) { return self.summary();})
            .def("put", +[](globimap_t &self, uint32_t x, uint32_t y) { self.put({x,y});})
            .def("get", +[](globimap_t &self, uint32_t x, uint32_t y) { return self.get({x,y});})
            .def("clear", +[](globimap_t &self){self.clear();})
            .def("map", +[](globimap_t &self, py::array_t<double> mat, int o0, int o1) {

                auto f = [&](int x, int y, double v){
                    if (v != 0 && v!= 1)
                        throw(std::runtime_error("data is not binary."));
                    if (v == 1)
                        self.put({static_cast<uint32_t>(o0+x),static_cast<uint32_t>(o1+y)});
                };

                map_matrix(mat, f);
            })
            .def("enforce", +[](globimap_t &self, py::array_t<double> mat, int o0, int o1){

                auto f = [&](int x, int y, double v){
                    if (v==0 && self.get({static_cast<uint32_t>(o0+x),static_cast<uint32_t>(o1+y)})) // this is a false positive
                        self.add_error({static_cast<uint32_t>(o0+x),static_cast<uint32_t>(o1+y)});
                };

                map_matrix(mat, f);
            })
            .def("rasterize", +[](globimap_t& self, size_t x, size_t y, size_t s0, size_t s1){
                std::vector<double>* v = &self.rasterize(x,y,s0,s1);

                std::unique_ptr<std::vector<double> > seq_ptr = std::make_unique<std::vector<double>>((*v)); // Unique ptr does not leak if for some reason py::capsule would throw.
                auto capsule = py::capsule(seq_ptr.get(), [](void *p) { std::unique_ptr<std::vector<double> >(reinterpret_cast<std::vector<double> *>(p)); });
                seq_ptr.release();

                return py::array({ s0, s1 }, v->data(), capsule);
            })

            .def("correct", +[](globimap_t& self, size_t x, size_t y, size_t s0, size_t s1) {
                std::vector<double>* v = &self.apply_correction(x, y, s0, s1);

                std::unique_ptr<std::vector<double> > seq_ptr = std::make_unique<std::vector<double>>((*v)); // Unique ptr does not leak if for some reason py::capsule would throw.
                auto capsule = py::capsule(seq_ptr.get(), [](void *p) { std::unique_ptr<std::vector<double> >(reinterpret_cast<std::vector<double> *>(p)); });
                seq_ptr.release();

                return py::array_t({ s0, s1 }, v->data(), capsule);
            });


    py::class_<globimap3d_t>(m, "globimap3d",py::module_local())
            .def(py::init<>())
            .def("configure", +[](globimap3d_t &self, size_t k, size_t m)  {  self.configure(k,m);})
            .def("summary", +[](globimap3d_t &self) { return self.summary();})
            .def("put", +[](globimap3d_t &self, uint64_t x, uint64_t y, uint64_t z) {
	    // arrange everything in first 128 bit of array
	    const uint64_t w = 0;
	    uint64_t a = (x << 32) | y, b = (z << 32) | w;
	    self.put({a,b});})
            .def("get", +[](globimap3d_t &self, uint64_t x, uint64_t y, uint64_t z) { 
	    const uint64_t w = 0;
	    uint64_t a = (x << 32) | y, b = (z << 32) | w;
	    return self.get({a,b});})
	    
            .def("clear", +[](globimap3d_t &self){self.clear();});
            /*  THESE ARE NOT YET 3D
		
.def("map", +[](globimap_t &self, py::array_t<double> mat, int o0, int o1) {

                auto f = [&](int x, int y, double v){
                    if (v != 0 && v!= 1)
                        throw(std::runtime_error("data is not binary."));
                    if (v == 1)
                        self.put({static_cast<uint32_t>(o0+x),static_cast<uint32_t>(o1+y)});
                };

                map_matrix(mat, f);
            })
            .def("enforce", +[](globimap_t &self, py::array_t<double> mat, int o0, int o1){

                auto f = [&](int x, int y, double v){
                    if (v==0 && self.get({static_cast<uint32_t>(o0+x),static_cast<uint32_t>(o1+y)})) // this is a false positive
                        self.add_error({static_cast<uint32_t>(o0+x),static_cast<uint32_t>(o1+y)});
                };

                map_matrix(mat, f);
            })
            .def("rasterize", +[](globimap_t& self, size_t x, size_t y, size_t s0, size_t s1){
                std::vector<double>* v = &self.rasterize(x,y,s0,s1);

                std::unique_ptr<std::vector<double> > seq_ptr = std::make_unique<std::vector<double>>((*v)); // Unique ptr does not leak if for some reason py::capsule would throw.
                auto capsule = py::capsule(seq_ptr.get(), [](void *p) { std::unique_ptr<std::vector<double> >(reinterpret_cast<std::vector<double> *>(p)); });
                seq_ptr.release();

                return py::array({ s0, s1 }, v->data(), capsule);
            })

            .def("correct", +[](globimap_t& self, size_t x, size_t y, size_t s0, size_t s1) {
                std::vector<double>* v = &self.apply_correction(x, y, s0, s1);

                std::unique_ptr<std::vector<double> > seq_ptr = std::make_unique<std::vector<double>>((*v)); // Unique ptr does not leak if for some reason py::capsule would throw.
                auto capsule = py::capsule(seq_ptr.get(), [](void *p) { std::unique_ptr<std::vector<double> >(reinterpret_cast<std::vector<double> *>(p)); });
                seq_ptr.release();

                return py::array_t({ s0, s1 }, v->data(), capsule);
            });*/



} // THE MODULE END
