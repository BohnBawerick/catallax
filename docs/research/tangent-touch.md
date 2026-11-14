# Tangent-touch prototype

Outcome: parked. Interior touches classify. Touches on a border do not.
When: 11 November 2026.

## What we were trying to achieve

A Boolean that carries on when a plane only touches a spline face. Today that cut refuses. The aim was to find every point where the plane touches without crossing and sort it into one of three kinds:

- an isolated touch, where the Boolean can continue;
- two branches crossing, which is a truly self-touching result and a legitimate final refusal;
- a degenerate contact, which gets a named refusal.

## Why it was proposed

It is step R2 ([archived plan, section 12.6](../archive/PLAN-2026-11.md#126-robustness-and-many-core-steps)). The reason on record: exact tangent contact is the one case the spline-cut proof cannot finish, and Catallax meets it more often than a floating-point kernel does, because designers place things flush and tangent on purpose. Floating-point checks cannot certify a double root. Exact elimination can.

## What was tried

For a plane against a rational patch: eliminate the two derivatives of the plane equation, isolate the roots, check each candidate point exactly, and classify it by the sign of an exact determinant. One measured round on 17 cases in five placements at 1, 2 and 4 workers, 255 requests in all, on the unchanged 20M allowance.

## Result

- Every copy gave the same answer at every worker count, with the same work.
- No request ran out of the 20M allowance.
- Interior touches classify correctly.

The gaps:

- a touch along a whole curve is not resolved;
- an ordinary touch on the border of the face is left undecided;
- the cylinder case was not measured, and cone, sphere and torus are not supported;
- nothing in the kernel or the app consumes the result, so it could not be shown to carry a Boolean through.

## Why it was parked

The promotion rule asked for the plane cases to be decided. Border touches and curve touches are not, and there was no consumer to prove the point of it. The default tangent refusal stays in force.

## What happens next

On 13 November the owner ruled for fewer refusals, with a bounded answer that is never shown as proved, and the bounded answer kind was added to the kernel. This roadmap makes tangent contact the first rule to use it, in milestone 1.1. The exact classifier here remains available for interior touches.

## Sources

- [prototypes/tangent-touch.md](../prototypes/tangent-touch.md): https://github.com/BohnBawerick/catallax/pull/403
- The bounded answer kind: https://github.com/BohnBawerick/catallax/pull/576
- [archived plan](../archive/PLAN-2026-11.md) section 12.6, row R2
