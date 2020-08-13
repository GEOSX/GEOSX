#include "BufferAllocator.hpp"
#include "common/DataTypes.hpp"

#ifdef GEOSX_USE_CHAI
namespace geosx
{
bool prefer_pinned_buffer = true;

void
setPreferPinned( bool p )
{
  prefer_pinned_buffer = p;
}

bool
getPreferPinned()
{
  return prefer_pinned_buffer;
}

}  // namespace geosx

#endif
