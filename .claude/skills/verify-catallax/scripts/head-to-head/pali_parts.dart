// Authored benchmark driver. Compile against the pinned, unmodified Pali Ijo package.
import 'dart:convert';
import 'dart:io';
import 'package:pali_ijo/pali_ijo.dart';
import 'package:vector_math/vector_math_64.dart';

Vector3 vector(dynamic value) => Vector3(
    (value[0] as num).toDouble(), (value[1] as num).toDouble(), (value[2] as num).toDouble());
double scalar(dynamic value) => (value as num).toDouble();

BrepSolid evaluate(List<dynamic> nodes) {
  final bodies = <BrepSolid>[];
  for (final node in nodes) {
    final a = node['args'] as List<dynamic>;
    final op = node['op'];
    late BrepSolid part;
    switch (op) {
      case 'box':
        final size = vector(a[0]);
        part = BrepXform.translate(BrepPrimitives.box(size.x, size.y, size.z), vector(a[1]));
      case 'cylinder':
        part = BrepPrimitives.cylinder(radius: scalar(a[0]), height: scalar(a[1]), origin: vector(a[2]), axis: vector(a[3]));
      case 'sphere':
        part = BrepPrimitives.sphere(radius: scalar(a[0]), center: vector(a[1]));
      case 'cone':
        part = BrepPrimitives.frustum(radius: scalar(a[0]), topRadius: scalar(a[1]), height: scalar(a[2]), origin: vector(a[3]), axis: vector(a[4]));
      case 'torus':
        part = BrepPrimitives.torus(center: Vector3.zero(), axis: Vector3(0, 0, 1), majorRadius: scalar(a[0]), minorRadius: scalar(a[1]));
      case 'rotate_translate':
        final q = a[1];
        final matrix = Quaternion(scalar(q[1]), scalar(q[2]), scalar(q[3]), scalar(q[0])).normalized().asRotationMatrix();
        part = BrepXform.affine(bodies[a[0]], BrepAffineTransform(
          matrix.transform(Vector3(1, 0, 0)), matrix.transform(Vector3(0, 1, 0)), matrix.transform(Vector3(0, 0, 1)), vector(a[2])));
      case 'union':
      case 'cut':
      case 'common':
        part = IjoOps.booleanSolid(bodies[a[0]], bodies[a[1]],
          op == 'union' ? IjoOps.unionOp : op == 'cut' ? IjoOps.subtractOp : IjoOps.intersectOp);
      case 'blend':
        final body = bodies[a[0]];
        final indices = <int>{};
        for (final point in a[1]) {
          final p = vector(point);
          final picks = <int>[];
          for (var i = 0; i < body.edges.length; i++) {
            if ((a[4] == null || (a[4] == 'circle' && body.edges[i].curve is BrepCircle)) &&
                body.edges[i].curve.closestPoint(p).distanceTo(p) <= 1e-7) picks.add(i);
          }
          if (picks.isEmpty || (a[4] == null && picks.length != 1)) throw StateError('Edge pick not unique: $point, $picks');
          indices.addAll(picks);
        }
        part = a[3] == true
          ? IjoOps.chamferSolid(body, indices.toList(), scalar(a[2]))
          : IjoOps.filletSolid(body, indices.toList(), radius: scalar(a[2]), inverse: false);
      case 'shell_top':
        final body = bodies[a[0]];
        final faces = BrepFaceSelect.extreme(Vector3(0, 0, 1)).resolve(body);
        part = IjoOps.hollowSolid(body, faces, Vector3(0, 0, 1), wallThicknessMm: scalar(a[1]));
      default:
        throw UnsupportedError('Operation $op is not expressible by this adapter');
    }
    bodies.add(part);
  }
  return bodies.last;
}

void main(List<String> args) {
  final request = jsonDecode(File(args[0]).readAsStringSync());
  final directory = Directory(args[1])..createSync(recursive: true);
  final variants = <String, dynamic>{};
  final result = <String, dynamic>{'status': 'completed', 'variants': variants,
    'edit_mode': 'complete target script rebuild, no incremental reuse',
    'edit_survival': {'status': 'not-run', 'reason': 'Pali Ijo ancestry IDs are process-local; no persistent rebuild naming contract'}};
  void checkpoint() => File('${directory.path}/phases.json').writeAsStringSync(jsonEncode(result));
  dynamic phase(Map<String, dynamic> row, String key, dynamic Function() action) {
    final clock = Stopwatch()..start();
    try {
      final value = action();
      row[key] = {'status': 'answered', 'wall_ms': clock.elapsedMicroseconds / 1000};
      return value;
    } catch (error) {
      row[key] = {'status': error is UnsupportedError || error is BrepUnsupported ? 'not-expressible' : 'command-failed',
        'wall_ms': null, 'refusal_upper_bound_ms': clock.elapsedMicroseconds / 1000,
        'detail': '$error'};
      return null;
    } finally {
      checkpoint();
    }
  }
  for (final label in ['construction', 'edit']) {
    final row = <String, dynamic>{};
    variants[label] = row;
    final part = phase(row, 'build', () => evaluate(request[label] as List<dynamic>)) as BrepSolid?;
    if (part == null) continue;
    final document = IjoDocument.solid(part, name: request['part']);
    phase(row, 'display', () {
      final mesh = IjoOps.tessellateClosed(part);
      File('${directory.path}/$label.stl').writeAsBytesSync(IjoStl.write(mesh));
    });
    final measured = phase(row, 'measurement', () {
      final p = IjoMeasure.ofSolid(part);
      return {'volume': p.volume, 'area': p.surfaceArea,
        'centroid': [p.centroid.x, p.centroid.y, p.centroid.z], 'is_exact': p.isExact};
    });
    if (measured != null) row['measurement']['claimed'] = measured;
    phase(row, 'step', () => File('${directory.path}/$label.step').writeAsBytesSync(IjoIo.encode(document, IjoFormat.step)));
    row['step']['expected_faces'] = part.faces.length;
    checkpoint();
  }
}
