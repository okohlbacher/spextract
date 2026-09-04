// GCC 13 has no std::ranges::to (arrives in GCC 14). Cluster-build shim for the mzPeak C++ library:
// the four `| std::ranges::to<C>()` sites are sed'ed to `| spx_to<C>()`; this header (injected with
// -include) supplies a minimal pipeable converter that recurses for nested containers.
#pragma once
#if __cplusplus >= 202002L   // inert for meson's C++17 sanity probe
#include <ranges>
#include <type_traits>
#include <utility>
template <class C> struct spx_to_tag {};
template <class C> constexpr spx_to_tag<C> spx_to() { return {}; }
template <class C, class R> C spx_to_impl(R&& r)
{
  C out;
  for (auto&& e : r)
  {
    using V = typename C::value_type;
    if constexpr (std::is_constructible_v<V, decltype(e)>) out.emplace_back(std::forward<decltype(e)>(e));
    else out.emplace_back(spx_to_impl<V>(std::forward<decltype(e)>(e)));
  }
  return out;
}
template <class R, class C>
  requires std::ranges::input_range<std::remove_cvref_t<R>>
C operator|(R&& r, spx_to_tag<C>) { return spx_to_impl<C>(std::forward<R>(r)); }
#endif
