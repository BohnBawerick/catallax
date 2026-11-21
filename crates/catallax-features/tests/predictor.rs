//! The predictor's small controls (Checkpoint 2 decisions, section 4.5): one
//! block and one fillet, with the block's position written in short or long
//! rationals. The tier changes with them. A fillet this small ends far under
//! the smaller cap, so it is never looked at, reads nothing, and is accepted
//! as it was.
use catallax_features::{
    Built, EvaluationOutcome, evaluate, evaluate_incremental_outcome,
    predictor::{self, Choice, Path, Tier},
};
use catallax_numbers::{Decider, OperationWork, Real, with_operation_budget};
use catallax_recipe::Recipe;
use catallax_topology::Loop;
use serde_json::json;

/// The application's allowance for one build.
const ALLOWANCE: u64 = 100_000_000;

/// A 20 by 20 by 10 block moved along x by `shift`, with one top edge rounded
/// at 1 mm.
fn rounded_block(shift: &str) -> Recipe {
    let value = json!({"format": "catallax-recipe/5", "name": "predictor control",
        "result": "round", "steps": [
            {"id": "b", "op": "block", "size": ["20", "20", "10"],
             "pose": {"rotation": ["1", "0", "0", "0"], "translation": [shift, "0", "0"]}},
            {"id": "round", "op": "fillet", "body": "b", "radius": "1", "edges": [
                {"between": [{"step": "b", "role": "end"},
                             {"step": "b", "role": {"side": {"segment": "x_max"}}}]}]}]});
    Recipe::parse(&value.to_string()).unwrap()
}

/// A shift of `bits` bits in the numerator and in the denominator: one over
/// three, written with a common factor the kernel does not divide out
/// because the two share none.
fn long_shift(bits: u32) -> String {
    let power = (0..bits).fold(String::from("1"), |digits, _| double(&digits));
    format!("{}/{}", add_one(&power), add_one(&add_one(&power)))
}

fn double(digits: &str) -> String {
    let mut out = Vec::new();
    let mut carry = 0;
    for digit in digits.bytes().rev() {
        let value = (digit - b'0') * 2 + carry;
        out.push(b'0' + value % 10);
        carry = value / 10;
    }
    if carry > 0 {
        out.push(b'0' + carry);
    }
    out.reverse();
    String::from_utf8(out).unwrap()
}

fn add_one(digits: &str) -> String {
    let mut out: Vec<u8> = digits.bytes().collect();
    for digit in out.iter_mut().rev() {
        if *digit == b'9' {
            *digit = b'0';
        } else {
            *digit += 1;
            return String::from_utf8(out).unwrap();
        }
    }
    out.insert(0, b'1');
    String::from_utf8(out).unwrap()
}

/// The build the application calls: the exact builders, then the rolling ball.
fn with_fitted_path(recipe: &Recipe) -> (EvaluationOutcome, OperationWork, Vec<Choice>) {
    let ((outcome, work), choices) = predictor::with_log(|| {
        with_operation_budget(ALLOWANCE, || {
            evaluate_incremental_outcome(recipe, None, &Default::default(), &mut Decider::default())
        })
    });
    assert!(!work.exhausted, "{work:?}");
    (outcome.unwrap(), work, choices)
}

/// The build with no fitted construction behind its fillet.
fn exact_only(recipe: &Recipe) -> (Built, OperationWork, Vec<Choice>) {
    let ((built, work), choices) = predictor::with_log(|| {
        with_operation_budget(ALLOWANCE, || evaluate(recipe, &mut Decider::default()))
    });
    assert!(!work.exhausted, "{work:?}");
    (built.unwrap(), work, choices)
}

/// The exact route under a cap is the exact route: the same body, the same
/// charged work and the same count. The attempt ends under the high tier's
/// cap, so no look is taken of it.
fn the_route_is_the_exact_one(shift: &str, tier: Tier, limbs: u64) {
    let recipe = rounded_block(shift);
    let (plain, plain_work, unread) = exact_only(&recipe);
    // With no fitted construction there is no cap and nothing is read.
    assert!(unread.is_empty(), "{unread:?}");
    let (outcome, work, choices) = with_fitted_path(&recipe);
    let [choice] = choices.as_slice() else {
        panic!("one fillet makes one attempt: {choices:?}");
    };
    assert_eq!((choice.step.as_str(), choice.path), ("round", Path::Exact));
    // Decisions 4.5: the quick look costs at most 1% of the operation's
    // work. This operation is a few hundred units, so it has none.
    assert!(choice.look.is_none(), "{choice:?}");
    assert!(choice.attempt_work < predictor::HIGH_COMPLEXITY_CAP);
    // Nothing was refunded or read, and the attempt's charges are the
    // request's.
    assert_eq!(
        (work.used, work.refunded, work.uncharged),
        (plain_work.used, 0, plain_work.uncharged)
    );
    assert_eq!(choice.attempt_work + block_work(&recipe), work.used);
    // What a look would have read of this block: two planes meet at every
    // edge, so D is 1, and the tier follows the length of its rationals.
    let look = look_at_the_block(&recipe);
    assert_eq!((look.degree, look.limbs), (1, limbs));
    assert_eq!(look.tier(), tier);
    assert!(look.read > 0 && look.read <= predictor::LOOK_READS);
    // What is accepted is decided as before: the one checker, on the same body.
    let EvaluationOutcome::Exact(checkpoint) = outcome else {
        panic!("the exact builders finish this fillet");
    };
    let built = checkpoint.built();
    assert!(built.model.check(built.result).is_empty());
    assert_eq!(
        built.model.faces_of_body(built.result).len(),
        plain.model.faces_of_body(plain.result).len()
    );
}

/// The recipe's block alone.
fn block_of(recipe: &Recipe) -> Built {
    let mut block: serde_json::Value = serde_json::from_str(&recipe.to_text()).unwrap();
    block["steps"].as_array_mut().unwrap().pop();
    block["result"] = json!("b");
    let block = Recipe::parse(&block.to_string()).unwrap();
    exact_only(&block).0
}

/// The look a fillet of every edge of the recipe's block would be given.
fn look_at_the_block(recipe: &Recipe) -> predictor::Look {
    let block = block_of(recipe);
    let mut edges = Vec::new();
    for face in block.model.faces_of_body(block.result) {
        for boundary in &block.model.face(face).loops {
            if let Loop::Coedges(uses) = block.model.loop_(*boundary) {
                for coedge in uses {
                    let edge = block.model.coedge(*coedge).edge;
                    if !edges.contains(&edge) {
                        edges.push(edge);
                    }
                }
            }
        }
    }
    assert_eq!(edges.len(), 12);
    predictor::of_fillet(&block.model, block.result, &edges, &Real::from_i64(1))
}

/// The charged work of the recipe's block alone.
fn block_work(recipe: &Recipe) -> u64 {
    let mut block: serde_json::Value = serde_json::from_str(&recipe.to_text()).unwrap();
    block["steps"].as_array_mut().unwrap().pop();
    block["result"] = json!("b");
    let block = Recipe::parse(&block.to_string()).unwrap();
    exact_only(&block).1.used
}

#[test]
fn short_rationals_are_the_low_tier_and_the_exact_route_is_unchanged() {
    the_route_is_the_exact_one("0", Tier::Low, 2);
    assert_eq!(Tier::Low.cap(), 20_000_000);
}

#[test]
fn long_rationals_are_the_high_tier_and_an_attempt_under_its_cap_stays_exact() {
    // 200 bits each way: four limbs for the numerator and four for the
    // denominator, so D x L^2 is 64, past the low tier's 32.
    the_route_is_the_exact_one(&long_shift(200), Tier::High, 8);
    assert_eq!(Tier::High.cap(), 2_000_000);
}
