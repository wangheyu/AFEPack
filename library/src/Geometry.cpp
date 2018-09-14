/**
 * @file   Geometry.cpp
 * @author Robert Lie
 * @date   Sun Apr 29 21:22:29 2007
 * 
 * @brief  
 * 
 * 
 */

#include "Geometry.h"

AFEPACK_OPEN_NAMESPACE

Geometry::Geometry() : ind(-1)
{}

Geometry::Geometry(const Geometry& g) :
  ind(g.ind),
     vtx(g.vtx),
     bnd(g.bnd)
{}

Geometry::~Geometry()
{}

Geometry& Geometry::operator=(const Geometry& g)
{
  ind = g.ind;
  vtx = g.vtx;
  bnd = g.bnd;
  return *this;
}

int Geometry::index() const
{
  return ind;
}

int& Geometry::index()
{
  return ind;
}

int Geometry::n_vertex() const
{
  return vtx.size();
}

int Geometry::n_boundary() const
{
  return bnd.size();
}

const std::vector<int>& Geometry::vertex() const
{
  return vtx;
}

std::vector<int>& Geometry::vertex()
{
  return vtx;
}

const std::vector<int>& Geometry::boundary() const
{
  return bnd;
}

std::vector<int>& Geometry::boundary()
{
  return bnd;
}

int Geometry::vertex(int i) const
{
  return vtx[i];
}

int& Geometry::vertex(int i)
{
  return vtx[i];
}

int Geometry::boundary(int i) const
{
  return bnd[i];
}

int& Geometry::boundary(int i)
{
  return bnd[i];
}

bool isSame(const Geometry& g0, const Geometry& g1)
{
  int i, j, k;
  j = g0.vtx.size();
  k = g1.vtx.size();
  if (j != k) return false;
  for (i = 0;i < k;i ++) {
    for (j = 0;j < k;j ++) {
      if (g0.vtx[i] == g1.vtx[j])
	break;
    }
    if (j == k) return false;
  }
  return true;
}

std::istream& operator>>(std::istream& is, Geometry& g)
{
  int i, j;
	
  is >> g.index();
  is >> i;
  g.vertex().resize(i);
  for (j = 0;j < i;j ++)
    is >> g.vertex()[j];
  is >> i;
  g.boundary().resize(i);
  for (j = 0;j < i;j ++)
    is >> g.boundary()[j];
  return is;
}

std::ostream& operator<<(std::ostream& os, const Geometry& g)
{
  int i, j;
	
  os << g.index() << "\n";
  i = g.vertex().size();
  os << i << "\t";
  for (j = 0;j < i;j ++)
    os << g.vertex()[j] << " ";
  os << "\n";
  i = g.boundary().size();
  os << i << "\t";
  for (j = 0;j < i;j ++)
    os << g.boundary()[j] << " ";
  os << "\n";
  return os;
}

///////////////////////////////////////////////////////////////////////////

GeometryBM::GeometryBM() : bm(0)
{}

GeometryBM::GeometryBM(const GeometryBM& g) :
  Geometry(g), bm(g.bm)
{}

GeometryBM::~GeometryBM()
{}

GeometryBM& GeometryBM::operator=(const GeometryBM& g)
{
  dynamic_cast<Geometry&>(*this) = g;
  bm = g.bm;
  return *this;
}

int GeometryBM::boundaryMark() const
{
  return bm;
}

int& GeometryBM::boundaryMark()
{
  return bm;
}

std::istream& operator>>(std::istream& is, GeometryBM& g)
{
  is >> dynamic_cast<Geometry&>(g);
  is >> g.bm;
  return is;
}

std::ostream& operator<<(std::ostream& os, const GeometryBM& g)
{
  os << dynamic_cast<const Geometry&>(g);
  os << g.bm << "\n";
  return os;
}

namespace __hidden_namespace {

#define MaxBits ( sizeof(unsigned) * CHAR_BIT )
#define IScale  \
  ((unsigned)((MaxBits <= 32) ? ~(0u) : (0xffffffff << (MaxBits - 32))))

  struct m_str {
    double x;
    double y;
    double z;
    unsigned int index[3];
    int table;
  };

  int cmp_indx(const void *a, const void *b) {    
    struct m_str *tmpa = (struct m_str*) a;
    struct m_str *tmpb = (struct m_str*) b;

    if (tmpa->index[0] > tmpb->index[0]) return  1;
    if (tmpa->index[0] < tmpb->index[0]) return -1;
  
    if (tmpa->index[1] > tmpb->index[1]) return  1;
    if (tmpa->index[1] < tmpb->index[1]) return -1;
  
    if (tmpa->index[2] > tmpb->index[2]) return  1;
    if (tmpa->index[2] < tmpb->index[2]) return -1;
  
    return 0;      
  }


  /*--------------------------------------------------------------------*/
  /* 2D Hilbert Space-filling curve                                     */
  /* 此处三维希尔伯特空间填充曲线的编码和解码的代码来源于互联网         */

  void hsfc2d(
              unsigned   coord[] , /* IN: Normalized integer coordinates */
              unsigned * nkey ,    /* IN: Word length of key */
              unsigned   key[] )   /* OUT: space-filling curve key */
  {
    static int init = 0 ;
    static unsigned char gray_inv[ 2 * 2 ] ;

    const unsigned NKey  = ( 2 < *nkey ) ? 2 : (*nkey) ;
    const unsigned NBits = ( MaxBits * NKey ) / 2 ;

    unsigned i ;
    unsigned char order[2+2] ;
    unsigned char reflect ;
  
    /* GRAY coding */

    if ( ! init ) {
      unsigned char gray[ 2 * 2 ] ;
      register unsigned k ;
      register unsigned j ;

      gray[0] = 0 ;
      for ( k = 1 ; k < sizeof(gray) ; k <<= 1 ) {
        for ( j = 0 ; j < k ; j++ ) gray[k+j] = k | gray[k-(j+1)] ;
      }
      for ( k = 0 ; k < sizeof(gray) ; k++ ) gray_inv[ gray[k] ] = k ;
      init = 1 ;
    }

    /* Zero out the key */

    for ( i = 0 ; i < NKey ; ++i ) key[i] = 0 ;

    order[0] = 0 ;
    order[1] = 1 ;
    reflect = ( 0 << 0 ) | ( 0 );

    for ( i = 1 ; i <= NBits ; i++ ) {
      const unsigned s = MaxBits - i ;
      const unsigned c = gray_inv[ reflect ^ (
                                              ( ( ( coord[0] >> s ) & 01 ) << order[0] ) |
                                              ( ( ( coord[1] >> s ) & 01 ) << order[1] ) ) ];

      const unsigned off   = 2 * i ;                   /* Bit offset */
      const unsigned which = off / MaxBits ;           /* Which word to update */
      const unsigned shift = MaxBits - off % MaxBits ; /* Which bits to update */

      /* Set the two bits */

      if ( shift == MaxBits ) { /* Word boundary */
        key[ which - 1 ] |= c ;
      }
      else {
        key[ which ] |= c << shift ;
      }

      /* Determine the recursive quadrant */

      switch( c ) {
      case 3:
        reflect ^= 03 ;
      case 0:
        order[2+0] = order[0] ;
        order[2+1] = order[1] ;
        order[0] = order[2+1] ;
        order[1] = order[2+0] ;
        break ;
      }
    }
  }

  /*--------------------------------------------------------------------*/
  /* 3D Hilbert Space-filling curve                                     */
  /* 此处三维希尔伯特空间填充曲线的编码和解码的代码来源于互联网         */

  void hsfc3d(unsigned   coord[] , /* IN: Normalized integer coordinates */
              unsigned * nkey ,    /* IN: Word length of 'key' */
              unsigned   key[] )   /* OUT: space-filling curve key */
  {
    static int init = 0 ;
    static unsigned char gray_inv[ 2*2*2 ] ;

    const unsigned NKey  = ( 3 < *nkey ) ? 3 : (*nkey) ;
    const unsigned NBits = ( MaxBits * NKey ) / 3 ;

    unsigned i ;
    unsigned char axis[3+3] ;
  
    /* GRAY coding */

    if ( ! init ) {
      unsigned char gray[ 2*2*2 ] ;
      register unsigned k ;
      register unsigned j ;

      gray[0] = 0 ;
      for ( k = 1 ; k < sizeof(gray) ; k <<= 1 ) {
        for ( j = 0 ; j < k ; j++ ) gray[k+j] = k | gray[k-(j+1)] ;
      }
      for ( k = 0 ; k < sizeof(gray) ; k++ ) gray_inv[ gray[k] ] = k ;
      init = 1 ;
    }

    /* Zero out the key */

    for ( i = 0 ; i < NKey ; ++i ) key[i] = 0 ;

    axis[0] = 0 << 1 ;
    axis[1] = 1 << 1 ;
    axis[2] = 2 << 1 ;

    for ( i = 1 ; i <= NBits ; i++ ) {
      const unsigned s = MaxBits - i ;
      const unsigned c = gray_inv[
                                  (((( coord[ axis[0] >> 1 ] >> s ) ^ axis[0] ) & 01 ) << 0 ) |
                                  (((( coord[ axis[1] >> 1 ] >> s ) ^ axis[1] ) & 01 ) << 1 ) |
                                  (((( coord[ axis[2] >> 1 ] >> s ) ^ axis[2] ) & 01 ) << 2 ) ];
      unsigned n ;

      /* Set the 3bits */

      for ( n = 0 ; n < 3 ; ++n ) {
        const unsigned bit   = 01 & ( c >> ( 2 - n ) );  /* Bit value  */
        const unsigned off   = 3 * i + n ;               /* Bit offset */
        const unsigned which = off / MaxBits ;           /* Which word */
        const unsigned shift = MaxBits - off % MaxBits ; /* Which bits */

        if ( MaxBits == shift ) { /* Word boundary */
          key[ which - 1 ] |= bit ;
        }
        else {
          key[ which ] |= bit << shift ;
        }
      }

      /* Determine the recursive quadrant */

      axis[3+0] = axis[0] ;
      axis[3+1] = axis[1] ;
      axis[3+2] = axis[2] ;

      switch( c ) {
      case 0:
        axis[0] = axis[3+2];
        axis[1] = axis[3+1];
        axis[2] = axis[3+0];
        break ;
      case 1:
        axis[0] = axis[3+0];
        axis[1] = axis[3+2];
        axis[2] = axis[3+1];
        break ;
      case 2:
        axis[0] = axis[3+0];
        axis[1] = axis[3+1];
        axis[2] = axis[3+2];
        break ;
      case 3:
        axis[0] = axis[3+2] ^ 01 ;
        axis[1] = axis[3+0] ^ 01 ;
        axis[2] = axis[3+1];
        break ;
      case 4:
        axis[0] = axis[3+2];
        axis[1] = axis[3+0] ^ 01 ;
        axis[2] = axis[3+1] ^ 01 ;
        break ;
      case 5:
        axis[0] = axis[3+0];
        axis[1] = axis[3+1];
        axis[2] = axis[3+2];
        break ;
      case 6:
        axis[0] = axis[3+0];
        axis[1] = axis[3+2] ^ 01 ;
        axis[2] = axis[3+1] ^ 01 ;
        break ;
      case 7:
        axis[0] = axis[3+2] ^ 01 ;
        axis[1] = axis[3+1];
        axis[2] = axis[3+0] ^ 01 ;
        break ;
      default:
        exit(-1);
      }
    }
  }

  /*--------------------------------------------------------------------*/


  void fhsfc2d(
                double     coord[] , /* IN: Normalized floating point coordinates */
                unsigned * nkey ,    /* IN: Word length of key */
                unsigned   key[] )   /* OUT: space-filling curve key */
  {
    const unsigned imax = IScale ;
    const double c0 = ( coord[0] <= 0 ) ? 0 : imax * coord[0] ;
    const double c1 = ( coord[1] <= 0 ) ? 0 : imax * coord[1] ;
    unsigned c[2] ;
    c[0] = (unsigned) (( imax < c0 ) ? imax : c0 );
    c[1] = (unsigned) (( imax < c1 ) ? imax : c1 );
    hsfc2d( c , nkey , key );
  }

  void fhsfc3d(double     coord[] , /* IN: Normalized floating point coordinates */
               unsigned * nkey ,    /* IN: Word length of key */
               unsigned   key[] )   /* OUT: space-filling curve key */
  {
    const double imax = ~(0u);
    unsigned c[3] ;
    c[0] = (unsigned) (coord[0] * imax) ;
    c[1] = (unsigned) (coord[1] * imax) ;
    c[2] = (unsigned) (coord[2] * imax) ;
    hsfc3d( c , nkey , key );
  }

  void hilbert(double *x, double *y, int *N, int *table)
  {
    unsigned int index[2];

    double extrx[2];
    double extry[2];
    struct m_str *s;
    
    int i;

    double temp[2]={0.0, 0.0};
    unsigned leng=2;

  
    s = (struct m_str*)malloc((*N)*sizeof(struct m_str ));
      
    extrx[0]=extrx[1]=x[0];
    extry[0]=extry[1]=y[0];
    table[0]=1;

    for(i=1;i<*N;i++) { 
      extrx[0]=std::min(extrx[0],x[i]);
      extrx[1]=std::max(extrx[1],x[i]);
      extry[0]=std::min(extry[0],y[i]);
      extry[1]=std::max(extry[1],y[i]);
      table[i]=i+1;
    }

    
    for(i=0;i<*N;i++) {
      temp[0]=(x[i]-extrx[0])/(extrx[1]-extrx[0]);
      temp[1]=(y[i]-extry[0])/(extry[1]-extry[0]);

      fhsfc2d(temp, &leng, index);  // 计算index， 张量
    
      s[i].x=x[i];
      s[i].y=y[i];
      s[i].index[0]=index[0];
      s[i].index[1]=index[1];
      s[i].table=table[i];

    }
  
    qsort(s,*N,sizeof(struct m_str), cmp_indx);

    for(i=0;i<*N;i++) {
      x[i]=s[i].x;
      y[i]=s[i].y;
      table[i]=s[i].table;
    }

    free(s);
  }

  void hilbert(double *x, double *y, double *z, int *N, int *table)
  {
    unsigned int index[3];

    double extrx[2];
    double extry[2];
    double extrz[2];
    struct m_str *s;
    
    int i;

    double temp[3]={0.0, 0.0, 0.0};
    unsigned leng=3;

  
    s = (struct m_str*)malloc((*N)*sizeof(struct m_str ));
      
    extrx[0]=extrx[1]=x[0];
    extry[0]=extry[1]=y[0];
    extrz[0]=extrz[1]=z[0];
    table[0]=1;

    for(i=1;i<*N;i++) { 
      extrx[0]=std::min(extrx[0],x[i]);
      extrx[1]=std::max(extrx[1],x[i]);
      extry[0]=std::min(extry[0],y[i]);
      extry[1]=std::max(extry[1],y[i]);
      extrz[0]=std::min(extrz[0],z[i]);
      extrz[1]=std::max(extrz[1],z[i]);
      table[i]=i+1;
    }

    
    for(i=0;i<*N;i++) {
      temp[0]=(x[i]-extrx[0])/(extrx[1]-extrx[0]);
      temp[1]=(y[i]-extry[0])/(extry[1]-extry[0]);
      temp[2]=(z[i]-extrz[0])/(extrz[1]-extrz[0]);

      fhsfc3d(temp,&leng,index);  // 计算index， 张量
    
      s[i].x=x[i];
      s[i].y=y[i];
      s[i].z=z[i];
      s[i].index[0]=index[0];
      s[i].index[1]=index[1];
      s[i].index[2]=index[2];
      s[i].table=table[i];

    }
  
    qsort(s,*N,sizeof(struct m_str), cmp_indx);

    for(i=0;i<*N;i++) {
    
      x[i]=s[i].x;
      y[i]=s[i].y;
      z[i]=s[i].z;
      table[i]=s[i].table;
    
    }

    free(s);
  }

}

void hsfc_renumerate(int n,
                     double * x,
                     double * y,
                     int * table) {
  __hidden_namespace::hilbert(x, y, &n, table);
  for (int i = 0;i < n;++ i) { table[i] -= 1; }
}

void hsfc_renumerate(int n,
                     double * x,
                     double * y,
                     double * z,
                     int * table) {
  __hidden_namespace::hilbert(x, y, z, &n, table);
  for (int i = 0;i < n;++ i) { table[i] -= 1; }
}

AFEPACK_CLOSE_NAMESPACE

/**
 * end of file
 * 
 */
