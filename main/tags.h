#ifndef tags_h
#define tags_h


// lhs tags: whether and how much of the lhs of the pde (ie the sparse system to be solved)
// must be updated each time step.
class tag_constant_lhs {};
class tag_part_constant_lhs {};
class tag_non_constant_lhs {};

class tag_expand_components {};
class tag_dont_expand_components {};

class tag_free_head {};
class tag_pivoting_head {};
class tag_clamped_head {};
class tag_swimming_head {};

class tag_position {};
class tag_angle {};

#endif


