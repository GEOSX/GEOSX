#include "KernelBase.hpp"
#include "kernelJITCompileCommands.hpp"

namespace geosx
{
namespace finiteElement
{

#if defined(GEOSX_USE_CUDA)
  constexpr bool compilerIsNVCC = true;
#else
  constexpr bool compilerIsNVCC = false;
#endif

jitti::CompilationInfo getCompilationInfo()
{
  jitti::CompilationInfo info;

  info.compileCommand = kernelJIT_COMPILE_COMMAND;
  info.compilerIsNVCC = compilerIsNVCC;
  info.linker = kernelJIT_LINKER;
  info.linkArgs = kernelJIT_LINK_ARGS;

  info.templateFunction = "geosx::finiteElement::fooBar";

  std::string const currentFile = __FILE__;
  info.headerFile = currentFile.substr( 0, currentFile.size() - ( sizeof( "kernelJIT.cpp" ) - 1 ) )
                                 + "KernelBase.hpp";

  return info;
}

} // namespace finiteElement
} // namespace geosx
