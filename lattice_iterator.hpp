#ifndef LATTICE_ITERATOR__HPP
#define LATTICE_ITERATOR__HPP

#include <iterator>


#define LATTICE_CONSTRUCTOR(NAME)                                       \
private:                                                                \
 T& lattice;                                                            \
 typename T::ZZ pos;                                                    \
public:                                                                 \
 NAME(T& lattice_, typename T::Site pos_)                               \
   : lattice(lattice_), pos(pos_) {};


#define LATTICE_CONSTRUCTOR_CONST(NAME)                                 \
private:                                                                \
 const T& lattice;                                                      \
 typename T::ZZ pos;                                                    \
public:                                                                 \
 NAME(const T& lattice_, typename T::Site pos_)                         \
   : lattice(lattice_), pos(pos_) {};



#define LATTICE_ITERATOR(NAME)                                          \
 NAME<T>& operator++()                                                  \
 {                                                                      \
   pos += 1;                                                            \
   return *this;                                                        \
 };                                                                     \
                                                                        \
 NAME<T> operator++(int)                                                \
 {                                                                      \
   NAME<T> orig =                                                       \
     NAME(lattice, pos);                                                \
   pos += 1;                                                            \
   return orig;                                                         \
 };                                                                     \
                                                                        \
 bool operator==(const NAME<T>& other) const                            \
 {                                                                      \
   return pos == other.pos;                                             \
 };                                                                     \
                                                                        \
 bool operator!=(const NAME<T>& other) const                            \
 {                                                                      \
   return pos != other.pos;                                             \
 };



template<typename T>
class ConstLatticeSiteIterator
  : public std::iterator<std::forward_iterator_tag, T>
{
  LATTICE_CONSTRUCTOR_CONST(ConstLatticeSiteIterator)
  LATTICE_ITERATOR(ConstLatticeSiteIterator)

  const typename T::Site& operator*() const
  { 
    return pos;
  };
  
  const typename T::Site* operator->() const
  {
    return &pos;
  };
};


template<typename T>
class LatticePhiIterator
  : public std::iterator<std::forward_iterator_tag, T>
{
  LATTICE_CONSTRUCTOR(LatticePhiIterator)
  LATTICE_ITERATOR(LatticePhiIterator)

  typename T::Phi& operator*()
  { 
    return lattice.get_phi(pos);
  };
};


template<typename T>
class ConstLatticePhiIterator
  : public std::iterator<std::forward_iterator_tag, T>
{
  LATTICE_CONSTRUCTOR_CONST(ConstLatticePhiIterator)
  LATTICE_ITERATOR(ConstLatticePhiIterator)

  const typename T::ConstPhi& operator*() const
  { 
    return lattice.get_phi(pos);
  };
};


template<typename T>
class ConstLatticePhiXYIterator
  : public std::iterator<std::forward_iterator_tag, T>
{
  LATTICE_CONSTRUCTOR_CONST(ConstLatticePhiXYIterator)
  LATTICE_ITERATOR(ConstLatticePhiXYIterator)

  typename T::ConstPhiXY operator*() const
  { 
    return lattice.get_phi_xy(pos);
  };
};





#endif
