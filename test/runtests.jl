using Polygons
using Base.Test

tests = ["definitions", "addition" ]
for t in tests
  test_file = "$(t).jl"
  print_with_color(:green, "* $test_file\n")
  include(test_file)
end
