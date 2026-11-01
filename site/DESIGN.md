---
name: Catallax + SolidDraft
description: A dark engineering page where exact geometry and its evidence do the persuading.
colors:
  ink: "#0b0c0d"
  graphite: "#131619"
  steel: "#262b30"
  slate: "#7b838b"
  fog: "#b4b0a6"
  bone: "#ede7da"
  oxide: "#ff5b1f"
  signal: "#a6e36b"
  datum: "#7cc7ff"
typography:
  display:
    fontFamily: "Bricolage Grotesque, Helvetica Neue, Arial, sans-serif"
    fontSize: "clamp(2.5rem, 5.6vw, 5.6rem)"
    fontWeight: 560
    lineHeight: 0.98
    letterSpacing: "-0.042em"
  hero:
    fontFamily: "Bricolage Grotesque, Helvetica Neue, Arial, sans-serif"
    fontSize: "clamp(3rem, 6.9vw, 8.4rem)"
    fontWeight: 560
    lineHeight: 0.92
    letterSpacing: "-0.05em"
  accent:
    fontFamily: "Instrument Serif, Georgia, serif"
    fontWeight: 400
    letterSpacing: "-0.02em"
  title:
    fontFamily: "Bricolage Grotesque, Helvetica Neue, Arial, sans-serif"
    fontSize: "clamp(1.8rem, 2.8vw, 2.6rem)"
    fontWeight: 560
    letterSpacing: "-0.04em"
  body:
    fontFamily: "Bricolage Grotesque, Helvetica Neue, Arial, sans-serif"
    fontSize: "1.0625rem"
    fontWeight: 400
    lineHeight: 1.55
  small:
    fontFamily: "Bricolage Grotesque, Helvetica Neue, Arial, sans-serif"
    fontSize: "0.9375rem"
    fontWeight: 400
    lineHeight: 1.5
  note:
    fontFamily: "Bricolage Grotesque, Helvetica Neue, Arial, sans-serif"
    fontSize: "0.875rem"
    fontWeight: 400
    lineHeight: 1.55
  label:
    fontFamily: "JetBrains Mono, ui-monospace, SF Mono, Menlo, monospace"
    fontSize: "0.72rem"
    fontWeight: 500
    lineHeight: 1.5
    letterSpacing: "0.08em"
  code:
    fontFamily: "JetBrains Mono, ui-monospace, SF Mono, Menlo, monospace"
    fontSize: "0.844rem"
    fontWeight: 400
    lineHeight: 1.72
rounded:
  xs: "4px"
  sm: "8px"
  md: "14px"
  lg: "18px"
  pill: "999px"
spacing:
  gutter: "clamp(20px, 4vw, 56px)"
  section: "clamp(110px, 14vw, 210px)"
  grid-gap: "clamp(32px, 5vw, 96px)"
components:
  button-primary:
    backgroundColor: "{colors.oxide}"
    textColor: "{colors.ink}"
    rounded: "{rounded.pill}"
    padding: "15px 24px"
  button-primary-hover:
    backgroundColor: "{colors.bone}"
    textColor: "{colors.ink}"
  button-ghost:
    textColor: "{colors.bone}"
    rounded: "{rounded.pill}"
    padding: "15px 24px"
  button-ghost-hover:
    backgroundColor: "{colors.bone}"
    textColor: "{colors.ink}"
  control-pill:
    textColor: "{colors.fog}"
    rounded: "{rounded.pill}"
    height: "36px"
  shape-active:
    backgroundColor: "{colors.bone}"
    textColor: "{colors.ink}"
    rounded: "{rounded.pill}"
  card:
    backgroundColor: "{colors.graphite}"
    textColor: "{colors.fog}"
    rounded: "{rounded.lg}"
    padding: "30px"
  code-window:
    backgroundColor: "{colors.graphite}"
    textColor: "{colors.fog}"
    rounded: "{rounded.md}"
    padding: "22px 24px 26px"
---

# Design System: Catallax + SolidDraft

Recorded from the previous site in October 2026, after its audit fixes. This is a reference for that design, including the retained comparison page and video renderers, not the published homepage layout. The bare-bones homepage and long-form pages are defined by `bare.css` and `read.css`. Both designs share `brand/tokens.css`, which the CSS and canvas renderer read. If this file and the tokens disagree, the tokens win.

## Overview

**Creative North Star: "The Inspection Bench"**

A part sits on a dark bench under one hot light, and everything around it is measured. The wireframe solids are that part. The labels, readouts and drawing marks are the instruments. Orange is the light: it lands on the single thing that matters in a view and nowhere else.

The page is dense with real material. It shows Rust from the repository, a real CLI drawing, committed benchmark data and a real build transcript, and it trusts a technical reader to look. The tone is a machined surface, hard-edged and exact, with one italic word as the only soft thing in a headline. Motion follows the token comment, "premium, rigid metal": fast out, no overshoot, nothing bouncing.

This name and the language around it were inferred from the shipped site. The owner has not chosen them.

**Key Characteristics:**
- One dark theme. There is no light mode.
- One accent color on near-neutral surfaces.
- Three typefaces, each with one job.
- Drawing marks as ornament: corner ticks, hairlines, a faint 96px grid.
- Everything readable and reachable before any JavaScript runs.

## Colors

The palette is a warm-to-cool gray stack with one hot accent and two quiet data colors.

### Primary
- **Oxide** (#ff5b1f): The hot light. It marks the primary action, the italic word in the hero, the progress bar, the anatomy rail, and the accent edges of the wireframes.

### Secondary
- **Signal** (#a6e36b): Numbers and strings in code, and one tool in the comparison bars.
- **Datum** (#7cc7ff): Type names in code, and one tool in the comparison bars.

### Neutral
- **Ink** (#0b0c0d): The page background.
- **Graphite** (#131619): Raised panels, which are cards and code windows.
- **Steel** (#262b30): Hairlines, borders and the drawing grid. It is never text.
- **Slate** (#7b838b): Labels, captions and code comments. It passes AA on both Ink and Graphite, at 5.09:1 and 4.72:1.
- **Fog** (#b4b0a6): Body copy.
- **Bone** (#ede7da): Headlines, emphasis and the filled hover state.

### Named Rules

**The One Hot Light Rule.** Oxide marks a single focal point per view. A second orange element in the same view should be a reason to remove one.

**The Check Both Surfaces Rule.** Any new secondary text color must clear 4.5:1 on Graphite as well as on Ink. Slate once passed on Ink alone and failed on every panel.

## Typography

**Display Font:** Bricolage Grotesque (with Helvetica Neue, Arial)
**Accent Font:** Instrument Serif Italic (with Georgia)
**Label/Mono Font:** JetBrains Mono (with ui-monospace, SF Mono, Menlo)

**Character:** A tight, slightly eccentric grotesque carries the argument, a single serif italic word adds the human note, and a mono face handles anything measured. All three are self-hosted and subset to Latin.

### Hierarchy
- **Hero** (560, `clamp(3rem, 6.9vw, 8.4rem)`, 0.92): The h1 only, set as three deliberate lines.
- **Display** (560, `clamp(2.5rem, 5.6vw, 5.6rem)`, 0.98): Section headings, balanced.
- **Title** (560, `clamp(1.8rem, 2.8vw, 2.6rem)`): Card titles.
- **Body** (400, 1.0625rem, 1.55): Prose, capped at about 58ch.
- **Small and Note** (400, 0.9375rem and 0.875rem): Lists, captions and footnotes.
- **Label** (500, 0.72rem, 0.08em tracking, uppercase by CSS): Section names, HUD readouts and table headers.
- **Code** (400, 0.844rem, 1.72): Rust excerpts. Long lines wrap.

### Named Rules

**The One Italic Word Rule.** A headline gets at most one Instrument Serif Italic word, and that word carries the idea. Never set a whole line in it.

**The Sentence Case Source Rule.** Write labels in sentence case and let CSS uppercase them. Hard-coded capitals read badly in screen readers and in copied text.

**The Mono Means Measured Rule.** JetBrains Mono is for code, data, readouts and labels. It is not a costume for prose.

## Layout

The page is one long column of full-width sections with a fluid gutter and very generous vertical padding. Most sections split into two columns on a 5:7 or 6:6 ratio, with the sticky or shorter side on the left.

The breakpoints follow the content. Navigation collapses to a Menu button at 1100px. Two-column layouts stack at 900px. The status pill and readouts drop at 560px. Tall headlines use `min-height`, never a fixed height, so the hero grows with its text and never clips it.

The kernel anatomy section is the one pinned moment. It pins for 300% of the viewport only when the viewport is at least 640px tall. Otherwise it stays a plain list, which also serves no-JS and high zoom.

**The Grow, Don't Clip Rule.** Copy stays in normal flow. Absolute positioning is for decoration only.

## Elevation & Depth

The system is flat. Depth comes from tonal steps, Ink to Graphite, and from 1px Steel hairlines. There are no shadows. The only layered effects are the fixed drawing grid behind the page and a faint static grain above it.

**The Flat Bench Rule.** Surfaces rest on the bench. A panel earns separation with a tone step and a hairline, never a shadow.

## Shapes

Two shape families carry the page. Controls are full pills, and panels have soft 14px to 18px corners. Small radii, 4px to 8px, are for the skip link, the logo and form controls. Wireframes and drawing ticks stay hard-edged, so softness belongs to the interface and precision to the geometry.

## Components

### Buttons
Confident pills with a fill that rises from the bottom edge.
- **Shape:** Full pill, 15px by 24px padding.
- **Primary:** Oxide fill with Ink text.
- **Ghost:** Steel hairline with Bone text.
- **Hover / Focus:** A Bone fill slides up over 480ms and the text turns Ink. Focus is a 2px Oxide ring offset by 3px.

### Control pills
Small mono pills for the Menu button, the pause control and the Ko-fi link. They use a Steel hairline and Fog text, are 36px tall, and grow to 44px under coarse pointers.

### Pause control
A control pill beside each moving element: the hero, the anatomy stage and both videos. All of them share one state, saved per browser, and they are hidden under reduced motion, where nothing moves on its own.

### Shape picker
A pill-shaped segmented control in the hero. The active shape fills with Bone, and `aria-pressed` carries the state.

### Cards / Containers
- **Corner Style:** Softly rounded, `{rounded.lg}`.
- **Background:** Graphite.
- **Shadow Strategy:** None. See Elevation & Depth.
- **Border:** 1px Steel, moving to Slate on hover.
- **Internal Padding:** 30px.

Under a fine pointer, cards tilt slightly and a soft glare follows the cursor.

### Code windows
Graphite panels with a Steel hairline and a mono file path as the title. There are no fake window controls. Excerpts come from the repository and must match it.

### Benchmark table
A native table with explicit roles, generated by `render/comparison.py`. Rows lay out as a grid, and each value carries a thin bar on top, in Oxide for Catallax and Slate for the rest. Under 900px the Tool header is visually hidden but stays available to screen readers.

### Navigation
A fixed bar that hides on scroll down and returns on scroll up or on keyboard focus. Links are uppercase mono with a thin Oxide underline that grows on hover. Under 1100px they move into a panel behind a Menu button, which Escape closes.

### Wireframe canvas
The signature element, drawn by `js/wire.js` in Canvas 2D. Edges fall into four styles: primary Bone, secondary Slate, accent Oxide, and dashed hidden lines. Depth fades the stroke. The solid builds edge by edge on arrival, redraws at no more than 30 fps while rotating, and stops drawing when nothing moves. It is always illustration, never kernel output.

## Do's and Don'ts

### Do:
- **Do** read every color from `brand/tokens.css`. The canvas uses the same tokens.
- **Do** animate opacity and transform only, using `--ease-out` or GSAP's `expo.out`.
- **Do** keep copy in normal flow, so it can grow under zoom, text spacing and small screens.
- **Do** put a pause control beside anything that moves on its own.
- **Do** pad SplitText line masks, so descenders and italic overhang stay inside the clip.
- **Do** link a claim to the evidence for it.

### Don't:
- **Don't** use `visibility` or `autoAlpha` for reveals. They remove content from the tab order and the accessibility tree.
- **Don't** run SplitText over an element containing a link. It strips the link's accessible name and breaks its underline.
- **Don't** hide text that has already painted in order to animate it in.
- **Don't** add shadows, fake window controls, pulsing status dots or marquees. They were removed on purpose.
- **Don't** number the sections.
- **Don't** add a kicker to a new heading that already names its section. The nine existing section labels stay because those headings are slogans and the labels say what each section is. This is a known exception and the design detector still reports it.
