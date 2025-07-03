#include <pybind11/pybind11.h>
#include<pybind11/numpy.h>
#include<pybind11/stl.h>

#include<iostream>
#include<array>

#include "murmur.hpp"
namespace py = pybind11;

/*
Bloom Filter Interface
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
		const char *p = static_cast<const char *>(data)+offset;
		for (int i=0; i < mat.itemsize(); i++)
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
          
	  k = _k; m=_m;
	  filter.resize(m);
	  for (size_t i=0; i < filter.size(); i++)
	      filter[i] = 0;
      }
    void insert(py::array mat, int axis=0){
	// axis not yet supported, maybe never
	std::vector<unsigned char> hash_data;
        
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

  void tobuffer(std::vector<char> &buf) const
    {
	buf.clear();
	char ch=0;
	  
	// push data
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

    void insert(const BloomFilterFacade &other){
	if (other.m != m || other.k != k)
	throw(std::runtime_error("Configuration Mismatch in union"));
	for (size_t i=0; i < other.filter.size(); i++)
	  if (other.filter[i]) filter[i] = true;
    }
    void intersect(const BloomFilterFacade &other){
	if (other.m != m || other.k != k)
	throw(std::runtime_error("Configuration Mismatch in intersect"));
	for (size_t i=0; i < other.filter.size(); i++)
	   filter[i] = other.filter[i] && filter[i];
    }
    
    void frombuffer(const std::vector<char> &buf)
    {
	size_t k=0;
	for (size_t i=0; i < buf.size(); i++)
	{
	    char ch = buf[i];
	    for(size_t j=0; j < 8; j++)
	      if (k < m)
	        filter[k++] = ch & (1<<j);
	
	}
    }


};




PYBIND11_MODULE(bgdbloomfilter, m) {
    /**
     * NOT FOR PRODUCTIVITY - DEVELOPMENT ONLY
     * */

    m.doc() = R"pbdoc(
        Bloom Filter - TUM BGDM Version
        -----------------------

        .. currentmodule:: python_example

        .. autosummary::
           :toctree: _generate

	bloomfilter
    )pbdoc";
    /*Bloom Filter*/
    py::class_<BloomFilterFacade>(m, "bloomfilter",py::module_local())
            .def(py::init<>())
            .def("configure", +[](BloomFilterFacade &self, size_t k, size_t m)  {  self.configure(k,m);})
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
		 self.configure(k,m);
                 self.frombuffer(std::vector<char>(data, data + length));
		})
	    .def(py::pickle(
			    [](const BloomFilterFacade &p) { // __getstate__
			      std::vector<char> out;
			      p.tobuffer(out);
			      return py::make_tuple(p.k, p.m,py::bytes(out.data(),out.size()));
                              },
			    [](py::tuple t) { // __setstate__
			      if (t.size() != 3)
				throw std::runtime_error("Invalid State in __setstate__");
			      BloomFilterFacade bf;
			      bf.configure(t[0].cast<int>(),t[1].cast<int>());
			      py::buffer_info info(py::buffer(t[2]).request());
			      const char *data = reinterpret_cast<const char *>(info.ptr);
			      size_t length = static_cast<size_t>(info.size);
			      bf.frombuffer(std::vector<char>(data, data + length));
			      return bf;
			    }
		      ))
	     .def("union", +[](BloomFilterFacade &self, BloomFilterFacade &other){
		self.insert(other);
	     })
	     .def("intersect", +[](BloomFilterFacade &self, BloomFilterFacade &other){
		self.intersect(other);
	     })
	     .def("get_k",+[](BloomFilterFacade &self){return self.k;})
	     .def("get_m",+[](BloomFilterFacade &self){return self.m;})
	    
	    ;



} // THE MODULE END
