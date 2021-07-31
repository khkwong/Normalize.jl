abstract type NormalizeMethod end

struct BoxCox <: NormalizeMethod end
struct YeoJohnson <: NormalizeMethod end
