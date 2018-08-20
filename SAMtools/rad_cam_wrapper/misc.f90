! Miscellaneous dummy routines
!#include "fintrf.h"

subroutine task_abort()
    call mexErrMsgIdAndTxt("dummy_functions:task_abort", &
        "task_abort called")
end subroutine

subroutine task_bcast_float4(rank_from, buffer, length)
    implicit none
    integer rank_from
    real(4) buffer(*)
    integer length
    call mexErrMsgIdAndTxt("dummy_functions:task_bcast_float4", &
        "task_bcast_float4 called")
end subroutine
