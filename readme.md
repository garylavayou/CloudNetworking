# Cloud Networking
## Dimensioning Network Slicing by Adjusting Resource Price
The objective is to maximize the net social welfare, at the same time each network slice want to maximize their own net profit. 

**Note**
1. Compared with the shortest path, longer path will be sold at a chipper price to encourage the virtual slice operator to use it, since the slice owner will pay for more links (and to guarantee the latency constraint, the physical owner should allocate more bandwidth to the slice owner.)
2. The pricing should distinguish different service type (e.g. different latency).
Usually, stringent latency constraint corresponds to high user utility. 
For example, considering file downloading, watching online videos and playing online games. 
Playing online games has most stringent delay constraint, it has highest utility, while file downloading almost has no delay constraint, so it has lowest utility.
3. The cost function of underlying resources add a penalty. 

## Fast Slice Reconfiguration

Fast slice reconfiuration is performed at flow arrival/departure time scale. 
The traffic variation at this time scale will be relatively small, and thus we can maximize the slice's porfit while keep the reconfiguration cost low.

We propose two fast slice reconfiguration schemes:
1. *Reconfigure Slice with Constant VNF capacity (RSCV)*: can change flow route, bandwith, and VNF instance assignment;
2. *Reconfigure Slice with Elastic VNF capacity (RSEV)*: in addtion to the ability of RSCV, we can further adjust the VNF instance's capacity.


## Hibrid Slice Reconfiguration

At long time scale, the traffic variation might be huge, so that the traffic demand and available resource in the slice might be mismathced.
As a result, the slice cannot achieve its optimal profit, without reconfigure the resource allocation for the slice.
Therfore, we popose hybrid slice reconfiguration scheme that bring slice dimensioning with fast slice reconfiguration together.
Slice dimensioning is scheduled at longer time scale to address the resource mismatch issue.
It could be scheduled either periodically (e.g. time period, event period, etc.), or triggered by slice condition, such as resource utilization, user data rate, profit level, etc.
In the hybrid scheme, slice dimensiong has been extended with consideration on reconfiguration cost.
So like the fast reconfigration schemes, the reconfiuration cost could also be reduced at a large extent.

# Configuration  
## Debug Control
* `DEBUG`: The global variable is used to enable debug code (See also [getstructfields](E:/workspace/Matlab/Projects/Language/getstructfields.m)). To define debug code, use following syntax:
    ```matlab {.line-numbers}
    global DEBUG;
    if ~isempty(DEBUG) && DEBUG
        debug_code statements;
    end
    ``` 
    Since `DEBUG` has global effect, to temporarily modify the debug behavior in the code, we can use the following code:
    ```matlab {.line-numbers}
    global DEBUG;
    old_debug = DEBUG;  % DEBUG might be empty.
    DEBUG = ture;       % to disable debug, set to false
    ...
    DEBUG = old_debug;  % recover DEBUG setting
    ```
* `INFO`: The global variable is used to control output of normal information. The usage is similar to `DEBUG`.
* `TRACE`: The global variable is used to enable conditional breakpoint.
